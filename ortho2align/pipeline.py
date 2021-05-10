import json
import os
import random
import time
from collections import defaultdict
from functools import partial
from pathlib import Path
from resource import getrusage, RUSAGE_SELF, RUSAGE_CHILDREN
from subprocess import DEVNULL, run
from tempfile import TemporaryDirectory

from tqdm import tqdm

from .benchmark import calc_ortholog_metrics, trace_orthologs
from .fitting import HistogramFitter, KernelFitter, ranking
from .genomicranges import (BaseGenomicRangesList, FastaSeqFile,
                            GenomicRangesAlignment, GenomicRangesList,
                            SequencePath)
from .parallel import (ExceptionLogger, NonExceptionalProcessPool,
                       TimeoutProcessPool)
from .parsing import ParserException, parse_annotation
from .utils import simple_hist, slplot


class CmdProgressBar(tqdm):
    # A solution from https://github.com/tqdm/tqdm/issues/509#issuecomment-410537096 for correct nesting.
    def __new__(cls, *args, **kwargs):
        try:
            cls._instances = tqdm._instances
        except AttributeError:
            pass

        instance = super().__new__(cls, *args, **kwargs)

        tqdm._instances = cls._instances
        return instance

    def __init__(self, cmd_hints, desc, disable=False):
        super().__init__(cmd_hints,
                         bar_format="{desc}: {n_fmt}/{total_fmt}, {elapsed}{postfix}]",
                         desc=desc,
                         postfix=cmd_hints[0],
                         initial=1,
                         disable=disable)

    def update(self):
        self.postfix = self.iterable[self.n]
        super().update()


def bg_from_inter_ranges(genes_filename, name_regex, sample_size, output_filename, seed, silent=False):
    random.seed(seed)

    cmd_hints = ['parsing the annotation...',
                 'sampling intergenic regions...',
                 'writing results to file...',
                 'finished.']
    with CmdProgressBar(cmd_hints,
                        'Sampling background from intergenic regions',
                        disable=silent) as pbar:

        with open(genes_filename, 'r') as infile:
            genes = parse_annotation(infile, name_regex=name_regex)

        pbar.update()

        inter_genes = genes.inter_ranges()
        if len(inter_genes) < sample_size:
            raise ValueError(f'The number of observations ({sample_size}) '
                             'must be less than the number of intergenic '
                             f'regions ({len(inter_genes)}) derived from genes')
        samples = BaseGenomicRangesList(random.sample(inter_genes, k=sample_size))

        pbar.update()

        with open(output_filename, 'w') as outfile:
            samples.to_bed6(outfile)

        pbar.update()


def bg_from_shuffled_ranges(genes_filename, genome_filename, name_regex, sample_size, output_filename, seed):
    with open(genes_filename, 'r') as infile:
        genes = parse_annotation(infile,
                                 name_regex=name_regex,
                                 sequence_file_path=genome_filename)
    if len(genes) < sample_size:
        raise ValueError(f'The number of observations ({sample_size}) '
                         'must be less than the number of genes '
                         f'({len(genes)})')
    shuffled_genes = genes.shuffle_inside_chrom(seed=seed)
    sample_shuffled_genes = shuffled_genes.sample_granges(n=sample_size, seed=seed)

    with open(output_filename, 'w') as outfile:
        sample_shuffled_genes.to_bed6(outfile)


def _align_two_ranges_blast(query, subject, **kwargs):
    try:
        return query.align_blast(subject, **kwargs)
    except Exception as e:
        raise ExceptionLogger(e, {'query': query,
                                  'subject': subject})


def _align_range_seqfile(query, seqfile, **kwargs):
    try:
        return query.align_blast_seqfile(seqfile, **kwargs)
    except Exception as e:
        raise ExceptionLogger(e, {'query': query,
                                  'seqfile': seqfile})


def _estimate_bg_for_single_query_blast(query, bg_ranges, word_size, output_name, observations=1000, seed=0):
    scores = list()
    try:
        for bg_range in bg_ranges:
            alignment = _align_two_ranges_blast(query, bg_range, word_size=word_size)
            alignment_scores = [hsp.score for hsp in alignment.HSPs]
            scores += alignment_scores
        score_size = len(scores)
        if score_size > observations:
            random.seed(seed)
            scores = random.sample(scores, observations)
        with open(output_name, 'w') as outfile:
            json.dump(scores, outfile)
        return output_name, score_size
    except Exception as e:
        raise ExceptionLogger(e, query)


def _estimate_bg_for_single_query_seqfile(query, seqfile, word_size, output_name, observations=1000, seed=0):
    try:
        alignment = _align_range_seqfile(query,
                                         seqfile,
                                         word_size=word_size,
                                         max_target_seqs=observations)
        scores = [hsp.score for hsp in alignment.HSPs]
        score_size = len(scores)
        if score_size > observations:
            random.seed(seed)
            scores = random.sample(scores, observations)
        with open(output_name, 'w') as outfile:
            json.dump(scores, outfile)
        return score_size
    except Exception as e:
        raise ExceptionLogger(e, query)


def estimate_background(query_genes_filename,
                        bg_ranges_filename,
                        query_genome_filename,
                        subject_genome_filename,
                        outdir,
                        cores=1,
                        word_size=6,
                        observations=1000,
                        seed=0,
                        query_name_regex=None,
                        bg_name_regex=None,
                        silent=False,
                        stats_filename=None):
    cmd_hints = ['reading annotations...',
                 'getting sample sequences...',
                 'aligning samples...',
                 'finished.']

    with CmdProgressBar(cmd_hints,
                        'Estimating background',
                        disable=silent) as pbar:

        with open(query_genes_filename, 'r') as infile:
            query_genes = parse_annotation(infile,
                                           sequence_file_path=query_genome_filename,
                                           name_regex=query_name_regex)
        with open(bg_ranges_filename, 'r') as infile:
            bg_ranges = parse_annotation(infile,
                                         sequence_file_path=subject_genome_filename,
                                         name_regex=bg_name_regex)

        pbar.update()

        total_query_genes = len(query_genes)

        query_chromsizes = FastaSeqFile(query_genome_filename).chromsizes
        subject_chromsizes = FastaSeqFile(subject_genome_filename).chromsizes

        with TemporaryDirectory() as tempdirname:
            tempdir = SequencePath(tempdirname)

            query_genes.get_fasta(outfileprefix='query_',
                                  outdir=tempdir,
                                  chromsizes=query_chromsizes)
            bg_seqfile = bg_ranges.get_fasta(outfileprefix='bg_seqfile',
                                             outdir=tempdir,
                                             chromsizes=subject_chromsizes,
                                             mode='bulk')

            pbar.update()

            if not os.path.exists(outdir):
                os.mkdir(outdir)
            data = ((query,
                     bg_seqfile,
                     word_size,
                     os.path.join(outdir, f"{query.name}.json"),
                     observations,
                     seed)
                    for query in query_genes)
            try:
                with NonExceptionalProcessPool(cores, verbose=(not silent)) as p:
                    results, exceptions = p.starmap_async(_estimate_bg_for_single_query_seqfile, data)
                total_exceptions = len(exceptions)
                if total_exceptions > 0:
                    pbar.write(f'{len(exceptions)} exceptions occured.')
                    for exception in exceptions:
                        pbar.write(str(exception))
                with open(os.path.join(outdir, 'exceptions.bed'), 'w') as outfile:
                    for exception in exceptions:
                        outfile.write(exception.variable.to_bed6() + '\n')
            except Exception as e:
                raise ExceptionLogger(e, p, 'Consider allocating more RAM or CPU time.')
        pbar.update()
        stats_msg = "-----------------------\n" \
                    f"estimate_background stats:\n" \
                    f"Recieved {total_query_genes} transcripts.\n" \
                    f"Estimated background for {total_query_genes - total_exceptions} of them.\n" \
                    f"Caught {total_exceptions} exceptions.\n" \
                    f"Distribution of amount of found background HSPs:\n{slplot(results)}\n" \
                    f"Reported at most {observations} HSPs for each transcript.\n" \
                    "-----------------------"
        if isinstance(stats_filename, str):
            with open(stats_filename, 'w') as stats_file:
                pbar.write(stats_msg, file=stats_file)
        else:
            pbar.write(stats_msg)


def _anchor(query_genes, query_anchors, query_genome_filename,
            subject_anchors, subject_genome_filename, query_chromsizes,
            subject_chromsizes, ortho_map, neighbour_dist, merge_dist, flank_dist):
    query_anchors.relation_mapping(subject_anchors,
                                   ortho_map,
                                   'orthologs')
    query_genes.get_neighbours(query_anchors,
                               neighbour_dist,
                               'neighbours')
    query_unalignable = list()
    query_prepared = list()

    for query_gene in query_genes:
        if len(query_gene.relations['neighbours']) == 0:
            query_unalignable.append(query_gene)
            continue
        syntenies = [grange
                     for neighbour in query_gene.relations['neighbours']
                     for grange in neighbour.relations['orthologs']]
        synteny_list = GenomicRangesList(syntenies,
                                         sequence_file_path=subject_genome_filename)
        query_gene.relations['syntenies'] = synteny_list.close_merge(merge_dist).flank(flank_dist,
                                                                                       chromsizes=subject_chromsizes)
        query_gene.relations['neighbours'] = query_gene.relations['neighbours'].close_merge(float('inf'))
        neighbourhood = query_gene.relations['neighbours'][0]
        if query_gene.end > neighbourhood.end or query_gene.start < neighbourhood.start:
            query_gene.relations['syntenies'] = query_gene.relations['syntenies'].flank(neighbour_dist + query_gene.end - query_gene.start,
                                                                                        chromsizes=subject_chromsizes)
        query_gene.relations['neighbours'] = GenomicRangesList([query_gene],
                                                               sequence_file_path=query_genome_filename)
        query_gene.relations['neighbours'] = query_gene.relations['neighbours'].flank(flank_dist,
                                                                                      chromsizes=query_chromsizes)

        query_prepared.append(query_gene)

    query_prepared_genes = GenomicRangesList(query_prepared,
                                             sequence_file_path=query_genome_filename)
    query_unalignable_genes = GenomicRangesList(query_unalignable,
                                                sequence_file_path=query_genome_filename)
    return query_prepared_genes, query_unalignable_genes


def _lift(query_genes, query_genome_filename, subject_genome_filename, query_chromsizes,
          subject_chromsizes, liftover_chains,
          merge_dist, flank_dist, min_ratio=0.05):

    query_unalignable = list()
    query_prepared = list()

    with TemporaryDirectory() as tempdirname:
        temp_query_filename = os.path.join(tempdirname, 'granges.bed')
        temp_result_filename = os.path.join(tempdirname, 'result.bed')
        temp_unmapped_filename = os.path.join(tempdirname, 'unmapped.bed')
        with open(temp_query_filename, 'w') as outfile:
            query_genes.to_bed6(outfile)
        run(f"liftOver -minMatch={min_ratio} -multiple -noSerial {temp_query_filename} {liftover_chains} {temp_result_filename} {temp_unmapped_filename}",
            shell=True,
            stderr=DEVNULL)
        with open(temp_result_filename, 'r') as infile:
            lifted = parse_annotation(infile)

    for query_gene in query_genes:
        if lifted.name_mapping.get(query_gene.name) is None:
            query_unalignable.append(query_gene)
            continue
        syntenies = [grange
                     for grange in lifted.name_mapping.get(query_gene.name)]
        synteny_list = GenomicRangesList(syntenies,
                                         sequence_file_path=subject_genome_filename)
        query_gene.relations['syntenies'] = synteny_list.close_merge(merge_dist).flank(flank_dist,
                                                                                       chromsizes=subject_chromsizes)
        query_gene.relations['neighbours'] = GenomicRangesList([query_gene],
                                                               sequence_file_path=query_genome_filename)
        query_gene.relations['neighbours'] = query_gene.relations['neighbours'].flank(flank_dist,
                                                                                      chromsizes=query_chromsizes)
        query_prepared.append(query_gene)

    query_prepared_genes = GenomicRangesList(query_prepared,
                                             sequence_file_path=query_genome_filename)
    query_unalignable_genes = GenomicRangesList(query_unalignable,
                                                sequence_file_path=query_genome_filename)
    return query_prepared_genes, query_unalignable_genes


def _align_syntenies(grange, word_size, outdir):
    try:
        neighbourhood = grange.relations['neighbours'][0]
        alignments = list()
        for synteny in grange.relations['syntenies']:
            alignment = neighbourhood.align_blast(synteny, word_size=word_size)
            alignment.to_genomic()
            alignment = alignment.cut_coordinates(qleft=grange.start,
                                                  qright=grange.end)
            alignment.qrange = grange
            alignments.append(alignment.to_dict())

        with open(os.path.join(outdir, f"{grange.name}.json"), 'w') as outfile:
            json.dump(alignments, outfile)

        return len(alignments)

    except Exception as e:
        raise ExceptionLogger(e, grange)


def get_alignments(mode,
                   query_genes,
                   query_genome,
                   subject_genome,
                   outdir,
                   query_anchors=None,
                   subject_anchors=None,
                   ortho_map=None,
                   liftover_chains=None,
                   query_name_regex=None,
                   query_anchors_name_regex=None,
                   subject_anchors_name_regex=None,
                   cores=1,
                   min_ratio=0.05,
                   neighbour_dist=0,
                   merge_dist=0,
                   flank_dist=0,
                   word_size=6,
                   silent=False,
                   stats_filename=None):

    program_mode = mode
    query_genes_filename = query_genes
    query_genome_filename = query_genome
    query_anchors_filename = query_anchors
    subject_anchors_filename = subject_anchors
    subject_genome_filename = subject_genome
    ortho_map_filename = ortho_map
    liftover_chains_filename = liftover_chains

    cmd_hints = ['reading annotations...',
                 'mapping orthology and neighbour relations, composing syntenies...',
                 'getting sequences...',
                 'getting alignments...',
                 'saving alignments...',
                 'finished.']

    with CmdProgressBar(cmd_hints,
                        'Getting alignments',
                        disable=silent) as pbar:

        with open(query_genes_filename, 'r') as infile:
            query_genes = parse_annotation(infile,
                                           sequence_file_path=query_genome_filename,
                                           name_regex=query_name_regex)
        total_query_genes = len(query_genes)
        query_genome = FastaSeqFile(query_genome_filename)
        query_chromsizes = query_genome.chromsizes
        subject_genome = FastaSeqFile(subject_genome_filename)
        subject_chromsizes = subject_genome.chromsizes

        if program_mode == "anchor":
            if any([arg is None
                    for arg in (query_anchors_filename,
                                subject_anchors_filename,
                                ortho_map_filename)]):
                raise ValueError('Arguments query_anchors, subject_anchors, ortho_map cannot be None for mode anchor.')
            with open(query_anchors_filename, 'r') as infile:
                query_anchors = parse_annotation(infile,
                                                 sequence_file_path=query_genome_filename,
                                                 name_regex=query_anchors_name_regex)

            with open(subject_anchors_filename, 'r') as infile:
                subject_anchors = parse_annotation(infile,
                                                   sequence_file_path=subject_genome_filename,
                                                   name_regex=subject_anchors_name_regex)

            with open(ortho_map_filename, 'r') as mapfile:
                ortho_map = json.load(mapfile)

            pbar.update()

            process_args = [query_genes, query_anchors, query_genome_filename,
                            subject_anchors, subject_genome_filename, query_chromsizes,
                            subject_chromsizes, ortho_map, neighbour_dist, merge_dist, flank_dist]

            query_prepared_genes, query_unalignable_genes = _anchor(*process_args)

        elif program_mode == "lift":
            pbar.update()
            if liftover_chains_filename is None:
                raise ValueError('Argument liftover_chains cannot be None for mode lift.')
            process_args = [query_genes, query_genome_filename, subject_genome_filename,
                            query_chromsizes, subject_chromsizes, liftover_chains_filename,
                            merge_dist, flank_dist, min_ratio]
            query_prepared_genes, query_unalignable_genes = _lift(*process_args)
        else:
            raise ValueError("program_mode is not one of ['anchor', 'lift']")

        total_unalignable_genes = len(query_unalignable_genes)

        pbar.update()

        query_chromsizes = query_prepared_genes.sequence_file.chromsizes

        with TemporaryDirectory() as tempdirname:
            tempdir = SequencePath(tempdirname)
            for query_gene in query_prepared_genes:
                query_gene.relations['neighbours'].get_fasta(outfileprefix='neigh_',
                                                             outdir=tempdir,
                                                             chromsizes=query_chromsizes)
                query_gene.relations['syntenies'].get_fasta(outfileprefix='synt_',
                                                            outdir=tempdir,
                                                            chromsizes=subject_chromsizes)

            pbar.update()

            if not os.path.exists(outdir):
                os.mkdir(outdir)

            align_syntenies_word_size = partial(_align_syntenies, word_size=word_size, outdir=outdir)

            with NonExceptionalProcessPool(cores, verbose=not silent) as p:
                try:
                    results, exceptions = p.map_async(align_syntenies_word_size, query_prepared_genes)
                except Exception as e:
                    raise ExceptionLogger(e, p, 'Consider allocating more RAM or CPU time.')

        query_exception_ranges = list()
        total_exceptions = len(exceptions)
        if total_exceptions > 0:
            pbar.write(f'{len(exceptions)} exceptions occured:')
            for exception in exceptions:
                pbar.write(str(exception))
                query_exception_ranges.append(exception.variable)

        pbar.update()

        query_exception_list = BaseGenomicRangesList(query_exception_ranges)

        with open(os.path.join(outdir, "unalignable.bed"), 'w') as outfile:
            query_unalignable_genes.to_bed6(outfile)

        with open(os.path.join(outdir, "exceptions.bed"), 'w') as outfile:
            query_exception_list.to_bed6(outfile)

        pbar.update()
        stats_msg = "-----------------------\n" \
                    f"get_alignments stats:\n" \
                    f"Recieved {total_query_genes} transcripts.\n" \
                    f"Aligned {total_query_genes - total_unalignable_genes - total_exceptions} of them.\n" \
                    f"Estimated as unalignable {total_unalignable_genes} transcripts.\n" \
                    f"Caught {total_exceptions} exceptions.\n" \
                    f"Distribution of amount of alignments:\n{simple_hist(results)}\n" \
                    f"Reported all alignments for each transcript.\n"\
                    "-----------------------"
        if isinstance(stats_filename, str):
            with open(stats_filename, 'a') as stats_file:
                pbar.write(stats_msg, file=stats_file)
        else:
            pbar.write(stats_msg)


def _build(query_gene_alignments_filename, query_bg_filename, fitter, fdr, pval_threshold):
    with open(query_gene_alignments_filename, 'r') as infile:
        query_gene_alignments = [GenomicRangesAlignment.from_dict(dict_alignment)
                                 for dict_alignment in json.load(infile)]
    if os.path.exists(query_bg_filename):
        with open(query_bg_filename, 'r') as infile:
            query_bg = json.load(infile)
    else:
        query_bg = [0, 1]  # simple workaround for correct KDE.
    if len(set(query_bg)) <= 1:
        query_bg += [1]
    if query_bg[0] == 1:
        query_bg.append(0)
    try:
        background = fitter(query_bg)
        score_threshold = background.isf(pval_threshold)
        aligned = False
        query_orthologs = list()
        subject_orthologs = list()
        query_dropped = list()
        for alignment in query_gene_alignments:
            # Applies to FDR also to reduce memory consumption.
            alignment.filter_by_score(score_threshold)
            # Heavily optimized FDR_BH procedure.
            if fdr:
                pvals = background.sf([hsp.score
                                       for hsp in alignment.HSPs])
                total_hsps = len(alignment._all_HSPs)
                retain = pvals <= total_hsps * pval_threshold / ranking(pvals)
                del pvals
                alignment.filter_by_bool_array(retain)
            alignment.srange.name = alignment.qrange.name
            result = alignment.best_chain()
            try:
                record = result.to_bed12(mode='str')
                query_orthologs.append(record[0])
                subject_orthologs.append(record[1])
                aligned = True
            except ValueError:
                query_dropped = [alignment.qrange]
        if aligned:
            query_dropped = list()
        return query_orthologs, subject_orthologs, query_dropped
    except Exception as e:
        raise ExceptionLogger(e, query_gene_alignments[0].qrange)


def build_orthologs(alignments,
                    background,
                    fitting,
                    outdir,
                    threshold=0.05,
                    fdr=False,
                    cores=1,
                    timeout=None,
                    silent=False,
                    stats_filename=None):

    alignments_dir = Path(alignments)
    background_dir = Path(background)
    fitting_type = fitting
    pval_threshold = threshold

    cmd_hints = ['reading alignments...',
                 'getting background...',
                 'refining alignments and calling orthologs...',
                 'saving to files...',
                 'finished.']

    with CmdProgressBar(cmd_hints,
                        'Building orthologs',
                        disable=silent) as pbar:

        alignments_files = [os.path.join(alignments_dir, filename)
                            for filename in os.listdir(alignments_dir)
                            if filename.endswith('json')]
        total_alignments = len(alignments_files)
        bg_files = [os.path.join(background_dir, os.path.basename(alignments_file))
                    for alignments_file in alignments_files]
        total_bgs = len(bg_files)

        pbar.update()

        if fitting_type == 'hist':
            fitter = HistogramFitter
        elif fitting_type == 'kde':
            fitter = KernelFitter

        pbar.update()

        build_partial = partial(_build,
                                fitter=fitter,
                                fdr=fdr,
                                pval_threshold=pval_threshold)
        data_for_building = zip(alignments_files, bg_files)

        with TimeoutProcessPool(cores, verbose=not silent) as p:
            try:
                orthologs, exceptions = p.starmap(build_partial,
                                                  data_for_building,
                                                  timeout)
            except Exception as e:
                raise ExceptionLogger(e, p, 'Consider allocating more RAM or CPU time.')

        total_exceptions = len(exceptions)
        query_exception_ranges = list()
        if total_exceptions > 0:
            pbar.write(f'{len(exceptions)} exceptions occured:')
            for exception in exceptions:
                pbar.write(str(exception))
                query_exception_ranges.append(exception.variable)

        pbar.update()

        if len(orthologs) == 0:
            query_orthologs, subject_orthologs, query_dropped = [], [], []
        else:
            query_orthologs, subject_orthologs, query_dropped = zip(*orthologs)

        dist_found = [len(group) for group in query_orthologs]
        query_orthologs = [ortholog
                           for group in query_orthologs
                           for ortholog in group]
        subject_orthologs = [ortholog
                             for group in subject_orthologs
                             for ortholog in group]
        query_dropped = BaseGenomicRangesList([grange
                                               for group in query_dropped
                                               for grange in group])
        total_dropped = len(query_dropped)
        query_exception_list = BaseGenomicRangesList(query_exception_ranges)

        if not os.path.exists(outdir):
            os.mkdir(outdir)
        query_output_filename = os.path.join(outdir, 'query_orthologs.bed')
        subject_output_filename = os.path.join(outdir, 'subject_orthologs.bed')
        query_dropped_filename = os.path.join(outdir, 'query_dropped.bed')
        query_exceptions_filename = os.path.join(outdir, 'query_exceptions.bed')

        with open(query_output_filename, 'w') as outfile:
            for record in query_orthologs:
                outfile.write(record + "\n")
        with open(subject_output_filename, 'w') as outfile:
            for record in subject_orthologs:
                outfile.write(record + "\n")
        with open(query_dropped_filename, 'w') as outfile:
            query_dropped.to_bed6(outfile)
        with open(query_exceptions_filename, 'w') as outfile:
            query_exception_list.to_bed6(outfile)

        pbar.update()
        stats_msg = "-----------------------\n" \
                    f"build_orthologs stats:\n" \
                    f"Recieved {total_alignments} transcript alignments and {total_bgs} backgrounds.\n" \
                    f"Built orthologs for {total_alignments - total_dropped - total_exceptions} of transcripts.\n" \
                    f"Found only unsignificant orthologs for {total_dropped} transcripts.\n" \
                    f"Caught {total_exceptions} exceptions.\n" \
                    f"Distribution of amount of orthologs:\n{simple_hist(dist_found)}\n" \
                    f"Reported all orthologs for each transcript.\n" \
                    "-----------------------"
        if isinstance(stats_filename, str):
            with open(stats_filename, 'a') as stats_file:
                pbar.write(stats_msg, file=stats_file)
        else:
            pbar.write(stats_msg)


def get_best_orthologs(query_orthologs,
                       subject_orthologs,
                       value,
                       function,
                       outfile_query,
                       outfile_subject,
                       outfile_map):

    query_filename = Path(query_orthologs)
    subject_filename = Path(subject_orthologs)
    value_name = value
    func_name = function
    query_output_filename = Path(outfile_query)
    subject_output_filename = Path(outfile_subject)
    outfile_map_filename = Path(outfile_map)

    extract = {'total_length': lambda line: (int(line.strip().split('\t')[2])
                                             - int(line.strip().split('\t')[1])),
               'block_count': lambda line: int(line.strip().split('\t')[9]),
               'block_length': lambda line: sum([int(i)
                                                 for i in line.strip().split('\t')[10].split(',')]),
               'weight': lambda line: float(line.strip().split('\t')[4]),
               'name': lambda line: line.strip().split('\t')[3]}

    func = {'min': min, 'max': max}

    name_pairs = list()

    with open(query_filename, 'r') as query_in, \
         open(subject_filename, 'r') as subject_in, \
         open(query_output_filename, 'w') as query_out, \
         open(subject_output_filename, 'w') as subject_out:
        purge = list()
        for query_line, subject_line in zip(query_in, subject_in):
            if len(purge) == 0:
                purge.append((query_line, subject_line))
            elif extract['name'](purge[-1][0]) == extract['name'](query_line):
                purge.append((query_line, subject_line))
            else:
                best_ortholog = func[func_name](purge,
                                                key=lambda i: extract[value_name](i[0]))
                query_out.write(best_ortholog[0])
                subject_out.write(best_ortholog[1])
                name_pairs.append([extract['name'](best_ortholog[0]),
                                   extract['name'](best_ortholog[1])])
                purge = list()
                purge.append((query_line, subject_line))

        if len(purge) > 0:
            best_ortholog = func[func_name](purge,
                                            key=lambda i: extract[value_name](i[0]))
            query_out.write(best_ortholog[0])
            subject_out.write(best_ortholog[1])
            name_pairs.append([extract['name'](best_ortholog[0]),
                               extract['name'](best_ortholog[1])])

    the_map = defaultdict(list)
    for query_name, subject_name in name_pairs:
        the_map[query_name].append(subject_name)
    with open(outfile_map_filename, 'w') as outfile:
        json.dump(the_map, outfile)


def annotate_orthologs(subject_orthologs,
                       subject_annotation,
                       output,
                       subject_name_regex=None,
                       stats_filename=None):

    subject_orthologs_filename = subject_orthologs
    subject_annotation_filename = subject_annotation
    output_filename = output

    try:
        with open(subject_orthologs_filename, 'r') as infile:
            subject_orthologs = parse_annotation(infile)
    except ParserException:
        with open(output_filename, 'w') as outfile:
            outfile.write('Query\tOrtholog\n')
        stats_msg = "-----------------------\n" \
                    f"annotate_orthologs stats:\n" \
                    f"Recieved {len(subject_orthologs)} orthologs."\
                    "-----------------------"
    else:
        with open(subject_annotation_filename, 'r') as infile:
            subject_annotation = parse_annotation(infile,
                                                  name_regex=subject_name_regex)

        subject_orthologs.get_neighbours(subject_annotation,
                                         relation='ortholog')

        name_map = {subject_ortholog.name: [grange.name
                                            for grange in subject_ortholog.relations['ortholog']]
                    for subject_ortholog in subject_orthologs}
        with open(output_filename, 'w') as outfile:
            outfile.write('Query\tOrtholog\n')
            dist_annot_amounts = list()
            for ortholog_name, annotation_names in name_map.items():
                if annotation_names:
                    outfile.write(f"{ortholog_name}\t{','.join(annotation_names)}\n")
                    dist_annot_amounts.append(len(annotation_names))
                else:
                    outfile.write(f"{ortholog_name}\tNotAnnotated\n")
                    dist_annot_amounts.append(0)

        stats_msg = "-----------------------\n" \
                    f"annotate_orthologs stats:\n" \
                    f"Recieved {len(subject_orthologs)} orthologs.\n" \
                    f"Distribution of amount of annotations:\n{simple_hist(dist_annot_amounts)}\n" \
                    f"Reported all annotations for each ortholog.\n" \
                    "-----------------------"
    if isinstance(stats_filename, str):
        with open(stats_filename, 'a') as stats_file:
            stats_file.write(stats_msg + '\n')
    else:
        print(stats_msg)


def benchmark_orthologs(query_genes,
                        query_name_regex,
                        found_query,
                        found_query_name_regex,
                        found_subject,
                        found_subject_name_regex,
                        found_query_map,
                        found_subject_map,
                        found_query_subject_map,
                        real_subject,
                        real_subject_name_regex,
                        real_map,
                        outfile,
                        tp_mode):

    query_genes_filename = query_genes
    found_query_filename = found_query
    found_subject_filename = found_subject
    found_query_map_filename = found_query_map
    found_subject_map_filename = found_subject_map
    found_query_subject_map_filename = found_query_subject_map
    real_subject_filename = real_subject
    real_map_filename = real_map
    out_filename = outfile

    with open(query_genes_filename, 'r') as infile:
        query_genes = parse_annotation(infile,
                                       name_regex=query_name_regex)

    with open(found_query_filename, 'r') as infile:
        found_query_genes = parse_annotation(infile,
                                             name_regex=found_query_name_regex)

    with open(found_subject_filename, 'r') as infile:
        found_subject_genes = parse_annotation(infile,
                                               name_regex=found_subject_name_regex)

    with open(real_subject_filename, 'r') as infile:
        real_subject_genes = parse_annotation(infile,
                                              name_regex=real_subject_name_regex)

    with open(found_query_map_filename, 'r') as infile:
        found_query_map = json.load(infile)

    with open(found_subject_map_filename, 'r') as infile:
        found_subject_map = json.load(infile)

    with open(found_query_subject_map_filename, 'r') as infile:
        found_query_subject_map = json.load(infile)

    with open(real_map_filename, 'r') as infile:
        real_map = json.load(infile)

    query_genes.relation_mapping(found_query_genes, found_query_map, 'found_query')
    query_genes.relation_mapping(found_subject_genes, found_subject_map, 'found_subject')
    query_genes.relation_mapping(real_subject_genes, real_map, 'real')

    for query_gene in query_genes:
        query_gene.relations['found_query'].relation_mapping(query_gene.relations['found_subject'],
                                                             found_query_subject_map,
                                                             'ortho_link')

    trace_orthologs(query_genes)
    metrics = calc_ortholog_metrics(query_genes, tp_mode=tp_mode)

    with open(out_filename, 'w') as outfile:
        json.dump(metrics, outfile)


def run_pipeline(query_genes,
                 query_genome,
                 subject_annotation,
                 subject_genome,
                 outdir,
                 query_anchors=None,
                 subject_anchors=None,
                 ortho_map=None,
                 liftover_chains=None,
                 query_anchors_name_regex=None,
                 subject_anchors_name_regex=None,
                 query_name_regex=None,
                 subject_name_regex=None,
                 sample_size=200,
                 observations=1000,
                 mode='lift',
                 min_ratio=0.01,
                 neighbour_dist=0,
                 merge_dist=0,
                 flank_dist=0,
                 fitting='kde',
                 threshold=0.05,
                 fdr=False,
                 timeout=None,
                 value='block_length',
                 function='max',
                 cores=1,
                 word_size=6,
                 seed=0,
                 silent=False,
                 annotate=False):

    start = time.time()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    bg_filename = os.path.join(outdir, 'shuffle_bg.bed')
    bg_outdir = os.path.join(outdir, 'bg_files')
    align_outdir = os.path.join(outdir, 'align_files')
    build_outdir = os.path.join(outdir, 'build_files')
    query_orthologs = os.path.join(build_outdir, 'query_orthologs.bed')
    subject_orthologs = os.path.join(build_outdir, 'subject_orthologs.bed')
    best_query_orthologs = os.path.join(outdir, 'best.query_orthologs.bed')
    best_subject_orthologs = os.path.join(outdir, 'best.subject_orthologs.bed')
    the_map = os.path.join(outdir, 'the_map.json')
    annotation_output = os.path.join(outdir, 'best.ortholog_annotation.tsv')
    stats_filename = os.path.join(outdir, 'stats.txt')

    bg_from_shuffled_ranges(genes_filename=subject_annotation,
                            genome_filename=subject_genome,
                            name_regex=subject_name_regex,
                            sample_size=sample_size,
                            output_filename=bg_filename,
                            seed=seed)
    estimate_background(query_genes_filename=query_genes,
                        bg_ranges_filename=bg_filename,
                        query_genome_filename=query_genome,
                        subject_genome_filename=subject_genome,
                        outdir=bg_outdir,
                        cores=cores,
                        query_name_regex=query_name_regex,
                        word_size=word_size,
                        observations=observations,
                        seed=seed,
                        silent=silent,
                        stats_filename=stats_filename)
    get_alignments(mode=mode,
                   query_genes=query_genes,
                   query_genome=query_genome,
                   subject_genome=subject_genome,
                   outdir=align_outdir,
                   query_anchors=query_anchors,
                   subject_anchors=subject_anchors,
                   ortho_map=ortho_map,
                   liftover_chains=liftover_chains,
                   query_name_regex=query_name_regex,
                   query_anchors_name_regex=query_anchors_name_regex,
                   subject_anchors_name_regex=subject_anchors_name_regex,
                   cores=cores,
                   min_ratio=min_ratio,
                   neighbour_dist=neighbour_dist,
                   merge_dist=merge_dist,
                   flank_dist=flank_dist,
                   word_size=word_size,
                   silent=silent,
                   stats_filename=stats_filename)
    build_orthologs(alignments=align_outdir,
                    background=bg_outdir,
                    fitting=fitting,
                    outdir=build_outdir,
                    threshold=threshold,
                    fdr=fdr,
                    cores=cores,
                    timeout=timeout,
                    silent=silent,
                    stats_filename=stats_filename)
    get_best_orthologs(query_orthologs=query_orthologs,
                       subject_orthologs=subject_orthologs,
                       value=value,
                       function=function,
                       outfile_query=best_query_orthologs,
                       outfile_subject=best_subject_orthologs,
                       outfile_map=the_map)
    if annotate:
        annotate_orthologs(subject_orthologs=best_subject_orthologs,
                           subject_annotation=subject_annotation,
                           output=annotation_output,
                           subject_name_regex=subject_name_regex,
                           stats_filename=stats_filename)
    end = time.time()
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(end - start))
    used_ram_mb = getrusage(RUSAGE_SELF).ru_maxrss / 1024 + getrusage(RUSAGE_CHILDREN).ru_maxrss / 1024
    with open(stats_filename, 'a') as stats_file:
        stats_file.write(f'Elapsed time: {elapsed_time}.\n')
        stats_file.write(f"Maximum RAM usage: {used_ram_mb // 1024:.0f} Gb {used_ram_mb % 1024:.0f} Mb.\n")

Search.setIndex({docnames:["index","modules","ortho2align"],envversion:{"sphinx.domains.c":1,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":1,"sphinx.domains.javascript":1,"sphinx.domains.math":2,"sphinx.domains.python":1,"sphinx.domains.rst":1,"sphinx.domains.std":1,sphinx:56},filenames:["index.rst","modules.rst","ortho2align.rst"],objects:{"":{ortho2align:[2,0,0,"-"]},"ortho2align.alignment":{Alignment:[2,1,1,""],HSP:[2,1,1,""],HSPVertex:[2,1,1,""],Transcript:[2,1,1,""],compare:[2,4,1,""],is_float:[2,4,1,""],numberize:[2,4,1,""],nxor_strands:[2,4,1,""]},"ortho2align.alignment.Alignment":{__eq__:[2,2,1,""],__init__:[2,2,1,""],_build_hsp_graph:[2,2,1,""],_find_all_transcripts:[2,2,1,""],_find_best_transcript:[2,2,1,""],_split_orientations:[2,2,1,""],best_transcript:[2,2,1,""],cut_coordinates:[2,2,1,""],filter_by_function_score:[2,2,1,""],filter_by_score:[2,2,1,""],from_dict:[2,2,1,""],from_file_blast:[2,2,1,""],get_all_transcripts:[2,2,1,""],get_best_transcripts:[2,2,1,""],plot_alignment:[2,2,1,""],replace_dict:[2,3,1,""],reset_filter_by_score:[2,2,1,""],set_zero:[2,2,1,""],to_dict:[2,2,1,""]},"ortho2align.alignment.HSP":{__eq__:[2,2,1,""],__init__:[2,2,1,""],_process_boundary:[2,2,1,""],copy:[2,2,1,""],distance:[2,2,1,""],from_dict:[2,2,1,""],orientation_dict:[2,3,1,""],precede:[2,2,1,""],to_dict:[2,2,1,""]},"ortho2align.alignment.HSPVertex":{__init__:[2,2,1,""],_relax:[2,2,1,""]},"ortho2align.alignment.Transcript":{HSPs:[2,3,1,""],__init__:[2,2,1,""],alignment:[2,3,1,""],from_dict:[2,2,1,""],plot_transcript:[2,2,1,""],score:[2,3,1,""],to_dict:[2,2,1,""]},"ortho2align.cli_scripts":{CLIVerbose:[2,1,1,""],align_syntenies:[2,4,1,""],align_two_ranges:[2,4,1,""],cache_orthodb_xrefs:[2,4,1,""],estimate_background:[2,4,1,""],get_alignments:[2,4,1,""],get_orthodb_map:[2,4,1,""],ortho2align:[2,4,1,""],refine_alignments:[2,4,1,""]},"ortho2align.cli_scripts.CLIVerbose":{__init__:[2,2,1,""]},"ortho2align.fitting":{AbstractFitter:[2,1,1,""],HistogramFitter:[2,1,1,""],KernelFitter:[2,1,1,""],inverse_approximator:[2,4,1,""]},"ortho2align.fitting.AbstractFitter":{__subclasshook__:[2,2,1,""],cdf:[2,2,1,""],estimator:[2,2,1,""],isf:[2,2,1,""],pdf:[2,2,1,""],ppf:[2,2,1,""],sf:[2,2,1,""]},"ortho2align.fitting.HistogramFitter":{__init__:[2,2,1,""],cdf:[2,2,1,""],data:[2,3,1,""],estimator:[2,3,1,""],isf:[2,2,1,""],pdf:[2,2,1,""],ppf:[2,2,1,""],sf:[2,2,1,""]},"ortho2align.fitting.KernelFitter":{__init__:[2,2,1,""],cdf:[2,2,1,""],data:[2,3,1,""],estimator:[2,3,1,""],isf:[2,2,1,""],pdf:[2,2,1,""],ppf:[2,2,1,""],sf:[2,2,1,""]},"ortho2align.genomicranges":{AlignedRangePair:[2,1,1,""],ChromosomeLocation:[2,1,1,""],ChromosomeNotFoundError:[2,5,1,""],EmptyGenomicRangesListError:[2,5,1,""],FastaSeqFile:[2,1,1,""],FilePathNotSpecifiedError:[2,5,1,""],GenomicCoordinatesError:[2,5,1,""],GenomicException:[2,5,1,""],GenomicRange:[2,1,1,""],GenomicRangesAlignment:[2,1,1,""],GenomicRangesList:[2,1,1,""],GenomicRangesTranscript:[2,1,1,""],InconsistentChromosomesError:[2,5,1,""],InconsistentGenomesError:[2,5,1,""],SequenceFileNotFoundError:[2,5,1,""],SequencePath:[2,1,1,""],align_with_relation_wrapper:[2,4,1,""],extract_taxid_mapping:[2,4,1,""],fasta_reformatter:[2,4,1,""]},"ortho2align.genomicranges.AlignedRangePair":{__eq__:[2,2,1,""],__init__:[2,2,1,""],__repr__:[2,2,1,""],__str__:[2,2,1,""],alignment:[2,3,1,""],from_dict:[2,2,1,""],query_grange:[2,3,1,""],subject_grange:[2,3,1,""],to_dict:[2,2,1,""]},"ortho2align.genomicranges.ChromosomeLocation":{size:[2,2,1,""],start:[2,2,1,""]},"ortho2align.genomicranges.FastaSeqFile":{__enter__:[2,2,1,""],__eq__:[2,2,1,""],__exit__:[2,2,1,""],__init__:[2,2,1,""],__repr__:[2,2,1,""],__str__:[2,2,1,""],_file_obj:[2,3,1,""],chromsizes:[2,3,1,""],from_dict:[2,2,1,""],get_fasta_by_coord:[2,2,1,""],line_length:[2,3,1,""],locate_coordinate:[2,2,1,""],sequence_file_path:[2,3,1,""],to_dict:[2,2,1,""]},"ortho2align.genomicranges.GenomicRange":{__eq__:[2,2,1,""],__init__:[2,2,1,""],__repr__:[2,2,1,""],__str__:[2,2,1,""],align_blast:[2,2,1,""],align_with_relation:[2,2,1,""],chrom:[2,3,1,""],distance:[2,2,1,""],end:[2,3,1,""],find_neighbours:[2,2,1,""],flank:[2,2,1,""],from_dict:[2,2,1,""],genome:[2,3,1,""],merge:[2,2,1,""],merge_relations:[2,2,1,""],name:[2,3,1,""],relations:[2,3,1,""],sequence_file_path:[2,3,1,""],sequence_header:[2,2,1,""],start:[2,3,1,""],strand:[2,3,1,""],to_dict:[2,2,1,""]},"ortho2align.genomicranges.GenomicRangesAlignment":{__init__:[2,2,1,""],from_dict:[2,2,1,""],from_file_blast:[2,2,1,""],to_dict:[2,2,1,""],to_genomic:[2,2,1,""],to_relative:[2,2,1,""]},"ortho2align.genomicranges.GenomicRangesList":{__eq__:[2,2,1,""],__init__:[2,2,1,""],__new__:[2,2,1,""],__repr__:[2,2,1,""],__str__:[2,2,1,""],_verbose_amount:[2,3,1,""],align_with_relation:[2,2,1,""],column_names:[2,3,1,""],comments:[2,3,1,""],dtypes:[2,3,1,""],end_types:[2,3,1,""],fileformats:[2,3,1,""],flank:[2,2,1,""],from_dict:[2,2,1,""],get_fasta:[2,2,1,""],get_neighbours:[2,2,1,""],inter_ranges:[2,2,1,""],key_func:[2,2,1,""],merge:[2,2,1,""],name_columns:[2,3,1,""],name_mapping:[2,3,1,""],name_patterns:[2,3,1,""],parse_annotation:[2,2,1,""],relation_mapping:[2,2,1,""],seps:[2,3,1,""],sequence_file:[2,3,1,""],sequence_file_path:[2,3,1,""],start_types:[2,3,1,""],to_dict:[2,2,1,""]},"ortho2align.genomicranges.GenomicRangesTranscript":{from_dict:[2,2,1,""],to_bed12:[2,2,1,""]},"ortho2align.genomicranges.SequencePath":{__eq__:[2,2,1,""],__fspath__:[2,2,1,""],__init__:[2,2,1,""],__repr__:[2,2,1,""],__str__:[2,2,1,""],check_correct:[2,2,1,""],exists:[2,2,1,""],open:[2,2,1,""],path:[2,3,1,""],to_dict:[2,2,1,""]},"ortho2align.orthodb":{FileOperationWrapper:[2,1,1,""],filter_odb_file:[2,4,1,""],filter_table:[2,4,1,""],load_odb_file:[2,4,1,""],load_table:[2,4,1,""],query_cached_odb_file:[2,4,1,""],query_table:[2,4,1,""],split_odb_file:[2,4,1,""],split_table:[2,4,1,""]},"ortho2align.orthodb.FileOperationWrapper":{__getitem__:[2,2,1,""],__init__:[2,2,1,""],__setitem__:[2,2,1,""],file_dict:[2,3,1,""],prefix:[2,3,1,""]},"ortho2align.parallel":{NonExceptionalProcessPool:[2,1,1,""],apply:[2,4,1,""],starapply:[2,4,1,""],starstarapply:[2,4,1,""]},"ortho2align.parallel.NonExceptionalProcessPool":{__enter__:[2,2,1,""],__exit__:[2,2,1,""],__init__:[2,2,1,""],_map_async:[2,2,1,""],executor:[2,3,1,""],map:[2,2,1,""],map_async:[2,2,1,""],max_workers:[2,3,1,""],shutdown:[2,2,1,""],starmap:[2,2,1,""],starmap_async:[2,2,1,""],suppress_exceptions:[2,3,1,""],verbose:[2,3,1,""]},ortho2align:{alignment:[2,0,0,"-"],cli_scripts:[2,0,0,"-"],fitting:[2,0,0,"-"],genomicranges:[2,0,0,"-"],orthodb:[2,0,0,"-"],parallel:[2,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","attribute","Python attribute"],"4":["py","function","Python function"],"5":["py","exception","Python exception"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:attribute","4":"py:function","5":"py:exception"},terms:{"abstract":2,"boolean":2,"byte":2,"case":2,"class":2,"default":2,"float":2,"function":2,"int":2,"new":2,"return":2,"static":2,"true":2,For:2,Has:2,NOT:2,One:2,The:2,There:2,These:2,Use:2,__enter__:2,__eq__:2,__exit__:2,__fspath__:2,__getitem__:2,__init__:2,__new__:2,__repr__:2,__setitem__:2,__str__:2,__subclasshook__:2,_all_hsp:2,_build_hsp_graph:2,_estim:2,_file_obj:2,_find_all_transcript:2,_find_best_transcript:2,_io:2,_map_async:2,_process_boundari:2,_relax:2,_split_orient:2,_verbose_amount:2,about:2,abstractfitt:2,accept:2,access:2,accordingli:2,add:2,added:2,addit:2,after:2,against:2,aim:2,algorithm:2,alia:2,align:1,align_blast:2,align_synteni:2,align_two_rang:2,align_with_rel:2,align_with_relation_wrapp:2,alignedrangepair:2,alignment_strand:2,all:2,allow:2,alreadi:2,ani:2,annot:2,api:2,appli:2,applier:2,approach:2,arg:2,arg_contain:2,argument:2,arithmet:2,aroudn:2,around:2,arrai:2,assign:2,associ:2,async:2,asynchron:2,attribut:2,auto:2,automat:2,avail:2,availab:2,averag:2,backward:2,bar:2,base:2,basic:2,bed12:2,bed3:2,bed6:2,befor:2,begin:2,beginn:2,being:2,best:2,best_transcript:2,between:2,biggest:2,binari:2,black:2,blast:2,blastn:2,block:2,blockcount:2,blocksiz:2,blockstart:2,blue:2,bool:2,both:2,boundari:2,build:2,built:2,bulk:2,cach:2,cache_orthodb_xref:2,cache_path:2,calcul:2,call:2,callabl:2,can:2,cast:2,cdf:2,cecond:2,certain:2,charact:2,cheaper:2,check:2,check_boundari:2,check_correct:2,choos:2,chose:2,chrom:2,chrom_start:2,chromosom:2,chromosome_nam:2,chromosomeloc:2,chromosomenotfounderror:2,chromosomseloc:2,chromsiz:2,classmethod:2,cli:2,cli_script:1,cliverbos:2,close:2,closer:2,closest:2,cls:2,code:2,coldtyp:2,collect:2,colnam:2,color:2,column:2,column_list:2,column_nam:2,combin:2,command:2,comment:2,compar:2,comparison:2,complement:2,complet:2,comput:2,concord:2,concurr:2,connect:2,consecut:2,consid:2,consist:2,contain:2,content:1,context:2,continu:2,conveni:2,coord:2,coordin:2,copi:2,core:2,coroutin:2,correct:2,correspond:2,cover:2,cpu_count:2,creat:2,criterion:2,cumul:2,current:2,cursor:2,custom:2,cut:2,cut_coordin:2,dash:2,data:2,databas:2,datafram:2,decreas:2,def:2,defalut:2,defin:2,definit:2,densiti:2,deriv:2,dfault:2,dict:2,dict_:2,dictionari:2,differ:2,direct:2,directori:2,distanc:2,distribut:2,divid:2,doe:2,doesn:2,down:2,dtype:2,dure:2,dynam:2,each:2,easier:2,edg:2,empir:2,empti:2,emptygenomicrangeslisterror:2,encod:2,end:2,end_typ:2,enter:2,epsilon:2,equal:2,estiam:2,estim:2,estimate_background:2,exc_tb:2,exc_typ:2,exc_val:2,exce:2,except:2,exclud:2,exclus:2,execut:2,executor:2,exist:2,exit:2,expect:2,expens:2,explan:2,explicitli:2,exponenti:2,extens:2,extract:2,extract_taxid_map:2,fals:2,fasta:2,fasta_reformatt:2,fastaseqfil:2,few:2,field:2,file:2,file_dict:2,file_object:2,file_path:2,file_suffix:2,fileformat:2,filenam:2,filenotspecifiederror:2,fileobj:2,fileoperationwrapp:2,filepathnotspecifiederror:2,filesystem:2,fill:2,filter:2,filter_by_function_scor:2,filter_by_scor:2,filter_odb_fil:2,filter_t:2,filtered_hsp:2,filtrat:2,find:2,find_neighbour:2,finish:2,first:2,fit:1,flag:2,flank:2,flank_dist:2,folder:2,follow:2,foo:2,foreign:2,form:2,format:2,found:2,fragment:2,from:2,from_dict:2,from_file_blast:2,full:2,func:2,futur:2,gap:2,gapextend:2,gapopen:2,gaussian_kd:2,gene:2,gene_id:2,geneid:2,gener:2,genom:2,genomic_rang:2,genomiccoordinateserror:2,genomicexcept:2,genomicrang:1,genomicrangelist:2,genomicrangesalign:2,genomicrangeslist:2,genomicrangestranscript:2,geq:2,get:2,get_align:2,get_all_transcript:2,get_best_transcript:2,get_fasta:2,get_fasta_by_coord:2,get_neighbour:2,get_orthodb_map:2,getter:2,gff:2,given:2,goe:2,grang:2,graph:2,greater:2,greedi:2,group:2,gtf:2,guarante:2,halt:2,handl:2,handler:2,have:2,header:2,helper:2,henc:2,here:2,high:2,histogram:2,histogramfitt:2,homolog:2,how:2,hsp:2,hsp_strand:2,hspvertex:2,ids:2,includ:2,inclus:2,inconsist:2,inconsistentchromosomeserror:2,inconsistentgenomeserror:2,increas:2,index:0,inessenti:2,inf:2,infer:2,infil:2,infin:2,inform:2,init:2,initi:2,inner:2,inplac:2,input:2,instal:2,instanc:2,integ:2,inter_rang:2,intersect:2,invalid:2,invers:2,inverse_approxim:2,invert:2,invok:2,is_float:2,isf:2,issubclass:2,item:2,itemrgb:2,iter:2,its:2,json:2,kde:2,kei:2,kept:2,kernelfitt:2,key_func:2,kwarg:2,lai:2,lambda:2,leav:2,left:2,length:2,leq:2,less:2,like:2,line:2,line_length:2,linewidth:2,link_color:2,list:2,list_of_chrom:2,load:2,load_odb_fil:2,load_tabl:2,locat:2,locate_coordin:2,lost:2,lowest:2,made:2,magic:2,mai:2,make:2,manag:2,mani:2,manner:2,manual:2,map:2,map_async:2,max_work:2,maximum:2,mean:2,memori:2,merg:2,merge_rel:2,messag:2,method:2,mimic:2,minim:2,minu:2,mode:2,modifi:2,modul:[0,1],monoton:2,more:2,multi:2,multipli:2,multiprocess:2,must:2,nai:2,name:2,name_column:2,name_map:2,name_pattern:2,ncbi:2,neg:2,neighbour:2,neq:2,next_vertic:2,non:2,none:2,nonexceptionalprocesspool:2,nonspecif:2,note:2,nucleotid:2,number:2,numer:2,nxor:2,nxor_strand:2,object:2,occur:2,odb_path:2,odb_prefix:2,odb_suffix:2,omit:2,one:2,one_chrom:2,one_genom:2,onli:2,open:2,oper:2,opposit:2,order:2,orient:2,orientation_dict:2,orthodb:1,ortholog:2,other:2,other_chrom:2,other_genom:2,other_list:2,otherwis:2,out:2,outfil:2,outfileprefix:2,outfmt:2,outprefix:2,output:2,overlap:2,own:2,packag:1,page:0,pair:2,panda:2,parallel:1,paramet:2,pars:2,parse_annot:2,parser:2,pass:2,path:2,pathlib:2,pattern:2,pdf:2,penalti:2,per:2,percent:2,perform:2,place:2,placement:2,plot:2,plot_align:2,plot_transcript:2,point:2,pointer:2,pool:2,posit:2,possibl:2,ppf:2,preced:2,preferr:2,prefix:2,present:2,previou:2,print:2,probabl:2,process:2,processpoolexecutor:2,produc:2,program:2,progress:2,properti:2,proportion:2,provid:2,providi:2,proxim:2,put:2,qchrom:2,qend:2,qleft:2,qlen:2,qlignment_strand:2,qname:2,qrang:2,qright:2,qstart:2,qstrand:2,queri:2,query_cached_odb_fil:2,query_grang:2,query_t:2,qzero:2,rais:2,rang:2,range_list_inst:2,read:2,read_csv:2,record:2,recov:2,red:2,reduc:2,refine_align:2,reformat:2,regardless:2,regex:2,region:2,rel:2,relat:2,relation_map:2,relax:2,replac:2,replace_dict:2,report:2,repres:2,represent:2,representaion:2,reset:2,reset_filter_by_scor:2,respect:2,restor:2,result:2,retriev:2,revers:2,right:2,rv_histogram:2,safe:2,same:2,scheme:2,schrom:2,scipi:2,score:2,search:[0,2],second:2,see:2,select:2,self:2,send:2,sep:2,separ:2,sequenc:2,sequence_fil:2,sequence_file_path:2,sequence_head:2,sequencefilenotfounderror:2,sequencepath:2,sequenti:2,set:2,set_scor:2,set_zero:2,setter:2,share:2,shortest:2,should:2,show:2,shut:2,shutdown:2,side:2,similar:2,simpl:2,singl:2,size:2,sleft:2,slen:2,smaller:2,smth1:2,smth2:2,sname:2,solid:2,some:2,sort:2,sortedcontain:2,sortedkeylist:2,sortedlist:2,sourc:2,speci:2,specif:2,specifi:2,split:2,split_odb_fil:2,split_tabl:2,srang:2,sright:2,sstart:2,sstrand:2,standalon:2,star:2,starappli:2,starmap:2,starmap_async:2,starstarappli:2,start:2,start_typ:2,statist:2,std:2,stdout:2,stop:2,store:2,str:2,strand:2,strategi:2,string:2,subjecct:2,subject:2,subject_grang:2,submodul:1,subtract:2,succeed:2,success:2,suffix:2,sum:2,support:2,suppress:2,suppress_except:2,surviv:2,sustain:2,symbol:2,sys:2,system:2,szero:2,tab:2,tabl:2,take:2,taken:2,target:2,task:2,taxid:2,ten:2,term:2,textiowrapp:2,than:2,thei:2,them:2,themselv:2,thi:2,thickend:2,thickstart:2,those:2,through:2,to_bed12:2,to_dict:2,to_genom:2,to_rel:2,togeth:2,total:2,toward:2,tqdm:2,transcript:2,transform:2,translat:2,tupl:2,turn:2,two:2,type:2,type_:2,undefin:2,under:2,uniqu:2,unless:2,until:2,usag:2,use:2,used:2,user:2,using:2,utf:2,util:2,valu:2,valueerror:2,verbos:2,veri:2,version:2,vertex:2,vertic:2,via:2,wai:2,wait:2,weight:2,well:2,were:2,when:2,where:2,whether:2,which:2,whitespac:2,within:2,won:2,wrapper:2,write:2,written:2,zero:2,zip:2},titles:["Welcome to ortho2align\u2019s documentation!","ortho2align","ortho2align package"],titleterms:{align:2,cli_script:2,content:2,document:0,fit:2,genomicrang:2,indic:0,modul:2,ortho2align:[0,1,2],orthodb:2,packag:2,parallel:2,submodul:2,tabl:0,welcom:0}})
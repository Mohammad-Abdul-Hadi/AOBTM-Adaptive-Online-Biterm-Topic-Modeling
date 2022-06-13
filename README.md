# Code of Adaptive Online Biterm Topic Model

  This framework takes reviews of different versions as input. To track the topic variations over versions, a novel method AOBTM (Adaptively Online Biterm Topic Modeling) is employed for generating version-sensitive topic distributions. The emerging topics are then identified based on the typical anomaly detection method. To make the topics comprehensible, the framework labels each topic with the most relevant phrases and sentences based on an effective ranking scheme considering semantic relevance. The prioritized topic labels are the app issues identified. Finally, it visualizes the variations of app issues along with versions, and highlights the emerging ones for better understanding. (70% done)
	
  The package contains Adaptive online algorithms for Biterm Topic Model (BTM): Adaptive Online BTM (AOBTM). AOBTM fits an individual BTM in a time slice by using the sufficient statistics as Dirichlet priors; it also takes window-size(w) into account to make the most recent version/time slice dependent on the previous w number of time-slices/versions.

More details can be referred to the following papers:

   > Xueqi Cheng, Xiaohui Yan, Yanyan Lan, and Jiafeng Guo. BTM: Topic Modeling over Short Texts. TKDE, 2014. 
   
   > Cuiyun Gao, Jichuan Zeng, Michael R. Lyu, and Irwin King. Online App Review Analysis for Identifying Emerging Issues. ICSE, 2018.
   
   > Mohammad A Hadi, Fatemeh H Fard, AOBTM: Adaptive Online Biterm Topic Model for Version Sensitive Short Text Analysis. ICSME 2020

## Usage ##

**Please use Python 2 for using this implemetation.**

The code includes a runnable example, you can run it by:

       $ script/runExample.sh

It trains BTM over the documents in *sample-data/0.txt, 1.txt, ...* and output the topics. The *n.txt* contains the training documents in time slice (supposed to be day) *n*, where each line represents one document with words separated by space as:
> word1 word2 word3 ....

(*Note: the sample data is only used for illustration of the usage of the code. It is not the data set used in the paper.*)

You can change the paths of data files and parameters in *script/runExample.sh* to run over your own data. 

Indeed, the *runExample.sh* processes the input documents in 4 steps.

**1. Index the words in the documents**   
   To simplify the main code, we provide a python script to map each word to a unique ID (starts from 0) in the documents. 

    $ python script/indexDocs.py <doc_dir> <dwid_dir> <voca_pt>
      doc_dir     input doc dir to be indexed, each file records docs in a day, while each line in a file is a doc with the format "word word ..."
      dwid_dir   output doc dir after indexing, each file records docs in a day, while each line is a doc with the format "wordId wordId ..."
      voca_pt   output vocabulary file, each line is a word with the format "wordId    word"

**2. Topic learning** 

       $ ./src/run obtm <K> <W> <alpha> <beta> <w> <n_iter> <docs_dir> <model_dir>
       or
       $ ./src/run ibtm <K> <W> <alpha> <beta> <n_iter> <docs_dir> <model_dir> <win> <n_rej>
       	K		int, number of topics
    	W		int, size of vocabulary
    	alpha		double, Symmetric Dirichlet prior of P(z), like 1
    	beta		double, Symmetric Dirichlet prior of P(w|z), like 0.01
		w		int, window-size (number of previous versions to depend on)
    	docs_dir    	string, path of training docs
    	model_dir	string, path of output directory
		win     	int, windows size of incremental Gibbs sampler
		n_rej   	int, rejuvenation sequence size of incremental Gibbs sampler

   The results will be written into the directory "model_dir":   
   - k20.pw_z: a K*M matrix for P(w|z), suppose K=20   
   - k20.pz:   a K*1 matrix for P(z), suppose K=20
   
**3. Inference topic proportions for documents, i.e., P(z|d)**     
   If you need to analysis the topic proportions of each documents, just run the following common to infer that using the model estimated.

    $ ./src/inf <type> <K> <day> <docs_dir> <model_dir>
      K	int, number of topics, like 20
      day   int, the nth day, like 0, 1, ..
      type	 string, 3 choices:sum_w, sum_b, mix. sum_b is used in our paper.
      docs_dir	string, path of training docs
      model_dir	string, output directory

   The result will be output to "model_dir":   
   - k20.day0.pz_d: a N*K matrix for P(z|d), suppose K=20 and day=0

**4. Results display**    
   Finally, we also provide a python script to illustrate the top words of the topics and their proportions in the collection. 

       $ python topicDisplay.py <model_dir> <K> <voca_pt>
	     model_dir    the output dir of BTM
	     K    the number of topics
	     voca_pt    the vocabulary file


## Citation ##
- The whole work is based on the following github repo:
> https://github.com/xiaohuiyan/OnlineBTM

If you have any questions, feel free to contact: [Mohammad Abdul Hadi](https://mohammad-abdul-hadi.github.io/).

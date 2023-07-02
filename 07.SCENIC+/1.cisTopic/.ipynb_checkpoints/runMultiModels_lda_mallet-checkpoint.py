import logging
import pickle
import sys
import argparse
import os
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from gensim import corpora, matutils

# arguments
def make_argument_parser():
	"""
	Creates an ArgumentParser to read the options for this script from
	sys.argv
	"""
	parser = argparse.ArgumentParser(description="Run topic models.")
	parser.add_argument('--mallet_path', '-m', type=str, required=False,
						default="/home/whe/Programs/package/Mallet-202108/bin/mallet",
						help='Path to mallet binary file.')
	parser.add_argument('--inputcisTopic_obj', '-i', type=str, required=True,
						help='Path to cisTopic object pickle file.')
	parser.add_argument('--output', '-o', type=str, required=True,
						help='Path to save final model list.')
	parser.add_argument('--n_topics', '-nt', type=str, required=True, nargs='+',
						help='Txt file containing selected topic id.')
	parser.add_argument('--n_cpu', '-c', type=int, required=True,
						help = 'Number of cores')
	parser.add_argument('--n_iter', '-it', type=int, required=False, default=150,
						help = 'Number of iterations')
	parser.add_argument('--alpha', '-a', type=int, required=False,  default=50,
						help='Alpha value')
	parser.add_argument('--alpha_by_topic', '-abt', type=str, required=False, default=True,
						help = 'Whether the alpha value should by divided by the number of topics')
	parser.add_argument('--eta', '-e', type=float, required=False, default=0.1,
						help='Eta value.')
	parser.add_argument('--eta_by_topic', '-ebt', type=str, required=False, default=False,
						help = 'Whether the eta value should by divided by the number of topics')
	parser.add_argument('--save_path', '-sp', type=str, required=False,
						default=None, help='Whether intermediate models should be saved')
	parser.add_argument('--seed', '-s', type=int, required=False,
						default=555, help='Seed for ensuring reproducibility')
	parser.add_argument('--temp_dir', '-td', type=str, required=False,
						default=None, help='Path to TMP directory')
	parser.add_argument('--init', '-init', default=True, action='store_true',
						help='whether initiate mallet')
	parser.add_argument('--no_init', '-no_init', dest='init', action='store_false',
						help='whether initiate mallet')
	parser.add_argument('--memory', '-mem', type=int, required=False,
						default=100, help='GB size for mallet memory')
	return parser


# multi topics
def run_cgs_multi_models_mallet(
	path_to_mallet_binary: str,
	cisTopic_obj: "cisTopicObject",
	n_topics: List[int],
	n_cpu: Optional[int] = 1,
	n_iter: Optional[int] = 150,
	random_state: Optional[int] = 555,
	alpha: Optional[float] = 50,
	alpha_by_topic: Optional[bool] = True,
	eta: Optional[float] = 0.1,
	eta_by_topic: Optional[bool] = False,
	top_topics_coh: Optional[int] = 5,
	tmp_path: Optional[str] = None,
	save_path: Optional[str] = None,
	reuse_corpus: Optional[bool] = False,
	init = True
):
	# Create cisTopic logger
	level = logging.INFO
	log_format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
	handlers = [logging.StreamHandler(stream=sys.stdout)]
	logging.basicConfig(level=level, format=log_format, handlers=handlers)
	log = logging.getLogger("cisTopic")

	# Init
	binary_matrix = cisTopic_obj.binary_matrix
	region_names = cisTopic_obj.region_names
	cell_names = cisTopic_obj.cell_names

	log.info(f"Formatting input to corpus")
	corpus = matutils.Sparse2Corpus(binary_matrix)
	names_dict = {x: str(x) for x in range(len(region_names))}
	id2word = corpora.Dictionary.from_corpus(corpus, names_dict)

	if init:
		init_model = LDAMallet(
			mallet_path=path_to_mallet_binary, id2word=id2word, num_topics=1, tmp_dir=tmp_path
		)
		init_model.convert_input(corpus, infer=False)
		del init_model

	# Run
	model_list = [run_cgs_model_mallet(
		path_to_mallet_binary, binary_matrix, corpus, id2word,
		n_topics=n,
		cell_names=cell_names,
		region_names=region_names,
		n_cpu=n_cpu,
		n_iter=n_iter,
		random_state=random_state,
		alpha=alpha,
		alpha_by_topic=alpha_by_topic,
		eta=eta,
		eta_by_topic=eta_by_topic,
		top_topics_coh=top_topics_coh,
		tmp_path=tmp_path,
		save_path=save_path,
		reuse_corpus=True
	) for n in n_topics]

	return model_list


# main	
def main():
	"""
	The main executable function
	"""
	
	parser = make_argument_parser()
	args = parser.parse_args()

	mallet_path = args.mallet_path
	print('Mallet path:', mallet_path)

	filename = args.inputcisTopic_obj
	infile = open(filename, 'rb')
	cisTopic_obj = pickle.load(infile)
	infile.close()
	print('Input cisTopic_object:', filename)
	
	output = args.output
	print('Output file:', output)
	
	n_topics = args.n_topics
	n_topics = list(map(int, n_topics[0].split(',')))
	print('Number of topics:', n_topics)
	
	alpha=args.alpha
	print('Alpha:', alpha)
	
	alpha_by_topic=args.alpha_by_topic
	print('Divide alpha by the number of topics:', alpha_by_topic)
	
	eta=args.eta
	print('Eta:', eta)
	
	eta_by_topic=args.eta_by_topic
	print('Divide eta by the number of topics:', eta_by_topic)
	
	n_iter=args.n_iter
	print('Number of iterations:', n_iter)
	
	n_cpu = args.n_cpu
	print('Number of cores:', n_cpu)

	save_path=args.save_path + "/"
	print('Path to save intermediate files:', save_path)
	if save_path == 'None':
		save_path = None
	
	random_state=args.seed
	print('Seed:', random_state)
	
	temp_dir=args.temp_dir + "/"
	print('Path to TMP dir:', temp_dir)

	init=args.init
	print('Whether initiate mallet:', str(init))

	mem=str(args.memory) + "G"
	print('Mallet memory:', mem)

	# Run models
	print('Running models')
	print('--------------')

	os.environ['MALLET_MEMORY'] = mem

	os.makedirs(save_path, exist_ok=True)
	os.makedirs(temp_dir, exist_ok=True)

	model_list=run_cgs_multi_models_mallet(
		mallet_path, cisTopic_obj,
		n_topics=n_topics,
		n_cpu=n_cpu,
		n_iter=n_iter,
		random_state=random_state,
		alpha=alpha,
		alpha_by_topic=alpha_by_topic,
		eta=eta,
		eta_by_topic=eta_by_topic,
		save_path=save_path,
		top_topics_coh=5,
		tmp_path=temp_dir,
		init=init
	)

	# Save
	with open(output, 'wb') as f:
		pickle.dump(model_list, f)

if __name__ == "__main__":
	main()

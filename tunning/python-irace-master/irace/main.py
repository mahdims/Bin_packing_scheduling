import argparse
import logging
import sys

def main(POP, CXPB, MUTPB, DATFILE):
	# just a test
	score = MUTPB*POP/100
	score = float(score)
	score = score - float(CXPB)
	if score < 0:
		score = 0

	# save the fo values in DATFILE
	with open(DATFILE, 'w') as f:
		f.write(str(score*100))

if __name__ == "__main__":
	# just check if args are ok
	with open('args.txt', 'w') as f:
		f.write(str(sys.argv))
	
	# loading example arguments
	ap = argparse.ArgumentParser(description='Feature Selection using GA with DecisionTreeClassifier')
	ap.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
	# 3 args to test values
	ap.add_argument('--pop', dest='pop', type=int, required=True, help='Population size')
	ap.add_argument('--cros', dest='cros', type=float, required=True, help='Crossover probability')
	ap.add_argument('--mut', dest='mut', type=float, required=True, help='Mutation probability')
	# 1 arg file name to save and load fo value
	ap.add_argument('--datfile', dest='datfile', type=str, required=True, help='File where it will be save the score (result)')

	args = ap.parse_args()
	logging.debug(args)
	# call main function passing args
	main(args.pop, args.cros, args.mut, args.datfile)
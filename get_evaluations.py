from evaluation.Program import Program
import sys

def main(seqfile):
    driver = Program(seqfile)
    driver.get_evaluation_result()


if __name__ == '__main__':
    main(sys.argv[1])
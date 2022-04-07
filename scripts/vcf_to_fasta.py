import argparse

def main():
    pass

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("contig")
    return parser

if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()

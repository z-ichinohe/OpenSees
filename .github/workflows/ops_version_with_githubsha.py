import argparse
parser: argparse.ArgumentParser = argparse.ArgumentParser()
parser.add_argument("github_sha", type=str)
args: argparse.Namespace = parser.parse_args()

with open("./SRC/OPS_Globals.h", "r", encoding="utf-8") as fp:
    lines: list[str] = fp.readlines()
    idx: int = [line.startswith("#define OPS_VERSION") for line in lines].index(True)
    line: str = lines[idx]
    line.replace('"\n', f'_{args.github_sha[0:7]}"\n')
    lines[idx] = line
with open("./SRC/OPS_Globals.h", "w", encoding="utf-8") as fp:
    fp.writelines(lines)
import argparse
import configparser
import os
import math
from decode import Decoder
from encode import Encoder
from utils.misc import check_integrity


def get_opts(init_config: str = "./configs/default.ini"):
    cp = configparser.ConfigParser()
    cp.read(init_config, encoding="UTF-8")

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--encode_file", default="./encode/50-SF.txt", type=str, help="input file path"
    )
    ap.add_argument(
        "--decode_file", default="./decode/50-SF.jpg", type=str, help="output file path"
    )
    ap.add_argument(
        "--origin_file",
        default="./origin/神奈川冲浪.jpg",
        type=str,
        help="original file path",
    )
    args = ap.parse_args()

    cp["DEFAULT"]["encode_file"] = args.encode_file
    cp["DEFAULT"]["decode_file"] = args.decode_file
    cp["DEFAULT"]["origin_file"] = args.origin_file

    return cp


def main():
    config = get_opts()
    origin_file = config["DEFAULT"]["origin_file"]
    encode_file = config["DEFAULT"]["encode_file"]
    decode_file = config["DEFAULT"]["decode_file"]
    need_encode = {"True": True, "False": False}[config["DEFAULT"]["need_encode"]]

    source_size = os.path.getsize(origin_file)
    chunk_size = 4
    chunk_num = math.ceil(source_size / chunk_size)

    if need_encode:
        print("Encoding...")
        Encoder(
            input_file=origin_file,
            output_file=encode_file,
            chunk_size=chunk_size,
            rs=2,
            max_homopolymer=3,
            gc=0.05,
            delta=0.001,
            c_dist=0.025,
            final=chunk_num * 5,
        ).encode()

    print("Decoding...")
    Decoder(
        input_file=encode_file,
        output_file=decode_file,
        chunk_num=chunk_num,
        header_size=4,
        rs=2,
        delta=0.001,
        c_dist=0.025,
        gc=0.05,
        max_homopolymer=3,
        max_hamming=0,
    ).decode()

    print(
        "Complete! With origin file: {}, encoded file: {}, recover file: {}".format(
            origin_file, encode_file, decode_file
        )
    )

    try:
        print("Checking integrity of output image...")
        check_integrity(origin_file, decode_file)
    except:
        print("Error occurred when trying to compare origin picture to output picture.")
        exit(1)


if __name__ == "__main__":
    main()

import logging
from typing import List, Tuple
from PIL import Image, ImageChops


def process_raw_input(input_file: str, chunk_size: int) -> Tuple[List[List[int]], int]:
    """
    read in file and transform to matrix
    :param
    input_file: file name; type: str
    chunk_size: number of information bytes per message; type: int

    :return
    data_array: data array; type: list.
    len_data: file size after zero padding; type: int.
    """

    with open(input_file, "rb") as file:
        data = file.read()

    # 原始文件二进制长度
    data_len = len(data)

    # 文件大小补齐为 chunk_size 倍数
    pad = -data_len % chunk_size
    if pad > 0:
        # Padded the file with zeros to have a round number of blocks of data
        data += b"\x00" * pad
        data_len = len(data)

    segments = data_len // chunk_size

    data_array: List[List[int]] = [None] * (segments)
    for num in range(segments):
        start = chunk_size * num
        end = chunk_size * (num + 1)
        chunk_binary = data[start:end]
        data_array[num] = [cb for cb in chunk_binary]

    return data_array, data_len


def check_integrity(origin_file: str, decode_file: str) -> None:
    origin_img = Image.open(origin_file)
    decode_img = Image.open(decode_file)

    diff = ImageChops.difference(origin_img, decode_img)
    logging.info("Checking integrity completed.")
    if diff.getbbox() is None:
        logging.info("No difference found between original image and decoded image!")
    else:
        logging.error("There are some differences between two images!")

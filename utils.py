from random import SystemRandom


def bin16_bin8():
    pass


def corrupt_chunk(chunk: list, n: int, k: int, t: int):
    cryptogen = SystemRandom()
    random_numbers = [cryptogen.randrange(n) for i in range(t)]
    random_positions = [cryptogen.randrange(len(chunk)) for i in range(k)]
    for _ in range(t):
        chunk[cryptogen.choice(random_positions)] = cryptogen.choice(random_numbers)
    return chunk


def char_to_int(msg_in: list):
    return [ord(char) for char in msg_in]


def int_to_char(msg_in: list):
    return ''.join(chr(char) for char in msg_in)


def chunk_message(msg_in: list, k: int):
    for i in range(0, len(msg_in), k):
        yield msg_in[i:i + k]

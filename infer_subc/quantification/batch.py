import os
import tempfile
from pandas import read_csv, concat
from pathlib import Path



# define function to append data to csv after each image is processed
def append_atomic_csv(csv_path, df):
    if isinstance(csv_path, str): csv_path = Path(csv_path)
    if df.empty:
        return
    if not os.path.exists(csv_path):
        df.to_csv(csv_path, index=False)
        return
    existing = read_csv(csv_path)
    merged = concat([existing, df], axis=0, ignore_index=True)
    merged = merged.drop_duplicates()  # final guard against repeat values
    fd, tmp_path = tempfile.mkstemp(prefix='append_', suffix='.csv')
    os.close(fd)
    merged.to_csv(tmp_path, index=False)
    os.replace(tmp_path, csv_path)


# define a function to check for existing data in csv files
def load_existing_keys_csv(csv_path, key_cols, chunksize=250_000):
    keys = set()
    if not os.path.exists(csv_path):
        return keys
    for chunk in read_csv(csv_path, usecols=key_cols, chunksize=chunksize):
        keys.update(map(tuple, chunk[key_cols].itertuples(index=False, name=None)))
    return keys
# ht1

This codebase implements a Python 3.6 transient detection algorithm for deployment on Loris NVIDIA Jetson Nano single-board computer. 

## Usage

1. Copy the script `ht1.py` from the `src/ht1` directory.
2. Run on HERMES SRA binary files with `python3 ht1.py SRAFILENAME`.

If the transient search returns a positive, a new file `SRAFILENAME_trigger.txt` is created in the `SRAFILENAME` parent directory. If no transient is found, no file is created.
The script returns error code 0 upon successful execution, error code 1 over wrong parameters, or error code 2 over invalid input file.

Search parameters can be controlled using a number of flags, see `python3 ht1.py --help` for more informations.

### Working with ht1 on ground

This package can also be installed for onground applications, see [documentation](docs/scripting_ht1.md) for more info.

## Testing
### With docker

The package comes with a minimal docker setup intended for testing in an replica environment of the NVIDIA Jetson Nano board.

1. Build the Docker container:
    ```bash
    docker-compose build
    ```
    Or, to build with no cache and logs, use `docker build --progress=plain --no-cache -t fm1trig-dev:latest . 2>&1`.
   
2. Run the container:
    ```bash
    docker-compose run --rm dev
    ```

3. Run the tests:
    ```bash
    python3 -m pytest tests/ -v
    ```

### Locally

Install the package with `pip install ".[dev]"` and run the tests with `python -m pytest tests`.
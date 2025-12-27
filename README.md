# FM1TRIG - SpIRIT GRB Trigger

Trigger algorithm for detecting gamma-ray bursts in SpIRIT satellite ratemeter data.

## Development Environment

This project uses Docker to replicate the target Ubuntu 18.04 / Python 3.6 environment.

### Prerequisites

- Docker
- Docker Compose

### Setup and Running

1. Build the Docker container:
    ```bash
    docker-compose build
    ```
    
    Or, to build with no cache and logs:
    
    ```bash
    docker build --progress=plain --no-cache -t fm1trig-dev:latest . 2>&1
    ```

2. Enter the container for interactive development:
    ```bash
    docker-compose run --rm dev
    ```

3. Optionally, install in editable mode.
   ```bash
   pip3 install -e ".[dev]"
   ```
   
### Running Tests

 ```bash
 python3 -m pytest tests/ -v
 ```

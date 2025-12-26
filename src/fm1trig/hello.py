#!/usr/bin/env python3
"""
Hello World script for fm1trig project.
Opens and reads a test SRA data file.
"""

import os

def hello_test():
    """Simple test function that returns True."""
    return True

def main():
    print("Hello from fm1trig!")

    # Path to a test data file
    test_file = "../../tests/data/20251018_054639/SRA/561DEC2C.raw"

    # Check if file exists
    if not os.path.exists(test_file):
        print(f"Error: Test file not found at {test_file}")
        return

    # Open and read the file
    print(f"\nOpening test file: {test_file}")
    with open(test_file, 'rb') as f:
        data = f.read()
        print(f"File size: {len(data)} bytes")
        print(f"First 64 bytes (hex): {data[:64].hex()}")

    print("\nSuccessfully read test file!")

if __name__ == "__main__":
    main()

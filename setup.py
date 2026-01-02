import pathlib

from setuptools import find_packages
from setuptools import setup

here = pathlib.Path(__file__).parent.resolve()
setup(
    name="ht1",
    version="0.1",
    author="Giuseppe Dilillo",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.6",
    extras_require={
        "dev": [
            "pytest",
            "pyinstrument",
        ],
    },
)

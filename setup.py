from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
setup(
    name="fm1trig",
    version="0.1",
    author="Giuseppe Dilillo",
    classifiers=[  # Optional
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires="~=3.6",
    extras_require={  # Optional
        "dev": ["pytest", "pyinstrument",],
    },
)
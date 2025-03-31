from setuptools import setup, find_packages

setup(
    name="mosplot",
    version="0.3.0",
    description="Mosfet plotting and simulation package.",
    author="Mohamed Watfa",
    author_email="medwatt@hotmail.com",
    url="https://github.com/medwatt/gmid",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)

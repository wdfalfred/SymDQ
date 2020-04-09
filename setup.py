import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="symdq", # Replace with your own username
    version="0.0.1",
    author="XIAO Yunchen",
    author_email="xiao.yc@foxmail.com",
    description="A package to do symbolic dual quaternion arithmetic.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ButteredCat/SymDQ",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
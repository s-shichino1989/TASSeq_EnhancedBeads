# distutils: extra_compile_args = ["-O3"]
# cython: language_level=3
__version__ = "1.10"
desc = 'part of BD Genomics Rhapsody Analysis pipeline version ' + __version__

def main():
    print(__version__, end="")

if __name__ == '__main__':
    main()

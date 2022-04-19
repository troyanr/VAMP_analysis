
from setuptools import setup, find_packages

setup(name='fastq_filters',
        version='1.0',
        description='Run Q score and mismatch filters on fastq files.',
        author='Aaron Maurais',
        url='https://github.com/ajmaurais/fastq_filters',
        classifiers=['Development Status :: 4 - Beta',
            'Intended Audience :: SCIENCE/RESEARCH',
            'Topic :: Software Development :: Build Tools',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            ],
        package_dir={'fastq_filters':'fastq_filters'},
        packages=find_packages(),
        python_requires='>=3.6.*',
        install_requires=['biopython==1.78'],
        entry_points={'console_scripts': ['q_score_filter=fastq_filters:run_q_score_filter',
                                          'mismatch_filter=fastq_filters:run_mismatch_filter',
                                          'print_alignment=fastq_filters:run_print_alignment']},
        )

# biosql_pyramid

Fork of https://github.com/chapmanb/biosqlweb which was build using pylons.
This project is ported to pyramid so you can run it locally without Google Apps Engine

# Setup

1. Clone repo

  ```
  $> git clone https://github.com/VDBWRAIR/biosql_pyramid
  $> cd biosql_pyramid
  ```

2. Get virtualenv setup(Optional, but recommended)

  ```
  $> virtualenv .
  $> . bin/activate
  ```

3. Install project and all dependencies

  ```
  $> python setup.py install
  ```

4. Initialize database

  ```
  $> ./initdb.py
  ```

5. Load some data from GenBank

  Create file called accession.lst with gi numbers or accessions(one per line)
  There is an example included with the following content

  ```
  KF824902
  KF824903
  JF504679
  ```

  Load it up

  ```
  $> ./loadseqs.py accessions.lst
  ``` 

6. Startup the server and have it open up browser

  ```
  pserve development.init -b
  ```
  
## Summary page

![Summary Page](/summary.png)

## Expanded Entry

![Expanded Entry](/expandedentry.png)

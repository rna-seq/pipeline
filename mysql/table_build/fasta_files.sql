CREATE TABLE IF NOT EXISTS fasta_files (
       species varchar(50) NOT NULL,
       genome varchar(100) NOT NULL,
       annotation varchar(100) NOT NULL,
       sidelength varchar(5) NOT NULL,
       filetype varchar(20) NOT NULL,
       location varchar(200) NOT NULL,
       index idx_genome (genome),
       index idx_annotation (annotation),
       index idx_filetype (filetype)
);

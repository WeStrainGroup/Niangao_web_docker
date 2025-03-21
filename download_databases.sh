#!/bin/bash

echo "Starting download_databases.sh..."

# 设置数据库下载目录
TAXONKIT_DIR="/opt/conda/envs/myenv/.taxonkit"
DB_DIR="/opt/conda/envs/myenv/db"
mkdir -p "$TAXONKIT_DIR"
mkdir -p "$DB_DIR"

# 下载并解压数据库
wget -P "$DB_DIR" https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_RefSeq_Fungi.tar.gz
tar -xzf "$DB_DIR/ITS_RefSeq_Fungi.tar.gz" -C "$DB_DIR"
rm "$DB_DIR/ITS_RefSeq_Fungi.tar.gz"

wget -P "$DB_DIR" https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_eukaryote_sequences.tar.gz
tar -xzf "$DB_DIR/ITS_eukaryote_sequences.tar.gz" -C "$DB_DIR"
rm "$DB_DIR/ITS_eukaryote_sequences.tar.gz"

wget -P "$DB_DIR" https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
tar -xzf "$DB_DIR/16S_ribosomal_RNA.tar.gz" -C "$DB_DIR"
rm "$DB_DIR/16S_ribosomal_RNA.tar.gz"

wget -P "$DB_DIR" https://ftp.ncbi.nlm.nih.gov/blast/db/28S_fungal_sequences.tar.gz
tar -xzf "$DB_DIR/28S_fungal_sequences.tar.gz" -C "$DB_DIR"
rm "$DB_DIR/28S_fungal_sequences.tar.gz"

wget -P "$DB_DIR" https://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences.tar.gz
tar -xzf "$DB_DIR/18S_fungal_sequences.tar.gz" -C "$DB_DIR"
rm "$DB_DIR/18S_fungal_sequences.tar.gz"

# 下载并解压 NCBI Taxonomy 数据库
wget -c -P "$DB_DIR" ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf "$DB_DIR/taxdump.tar.gz" -C "$DB_DIR"
rm "$DB_DIR/taxdump.tar.gz"

# 复制 TaxonKit 所需的文件
cp "$DB_DIR/names.dmp" "$DB_DIR/nodes.dmp" "$DB_DIR/delnodes.dmp" "$DB_DIR/merged.dmp" "$TAXONKIT_DIR"

echo "BLAST databases downloaded and decompressed to $DB_DIR"
echo "TaxonKit data copied to $TAXONKIT_DIR"
echo "download_databases.sh completed."
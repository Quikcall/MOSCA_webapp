CREATE DATABASE `myproject` DEFAULT CHARACTER SET latin1;

USE `myproject`;

DROP TABLE IF EXISTS `projects`;
DROP TABLE IF EXISTS `samples`;
DROP TABLE IF EXISTS `exe_projects`;
DROP TABLE IF EXISTS `users`;

CREATE TABLE `projects` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `project_name` varchar(50) DEFAULT NULL,
  `author` varchar(300) DEFAULT NULL,
  `description` text,
  `database_dir` varchar(300) DEFAULT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `threads` int(11) DEFAULT NULL,
  `sequencing` varchar(50) DEFAULT NULL,
  `quality_scores` varchar(50) DEFAULT NULL,
  `output_lvl` varchar(50) DEFAULT NULL,
  `data_type` varchar(50) DEFAULT NULL,
  `preprocessing` tinyint(1) DEFAULT NULL,
  `assembly` tinyint(1) DEFAULT NULL,
  `binning` tinyint(1) DEFAULT NULL,
  `assembler` varchar(50) DEFAULT NULL,
  `memory` varchar(20) DEFAULT NULL,
  `k_mer_sizes` varchar(300) DEFAULT NULL,
  `m` int(11) DEFAULT NULL,
  `alignment_method` varchar(50) DEFAULT NULL,
  `alignment_options` varchar(50) DEFAULT NULL,
  `train` varchar(300) DEFAULT NULL,
  `up_names_tax` varchar(2000) DEFAULT NULL,
  `up_sequences` varchar(500) DEFAULT NULL,
  `up_function` varchar(500) DEFAULT NULL,
  `up_miscellaneous` varchar(500) DEFAULT NULL,
  `up_interaction` varchar(500) DEFAULT NULL,
  `up_expression` varchar(500) DEFAULT NULL,
  `up_gene_ont` varchar(500) DEFAULT NULL,
  `up_chebi` varchar(500) DEFAULT NULL,
  `up_path_biot` varchar(500) DEFAULT NULL,
  `up_cell_loc` varchar(500) DEFAULT NULL,
  `up_ptm` varchar(500) DEFAULT NULL,
  `up_structure` varchar(500) DEFAULT NULL,
  `up_pubs` varchar(500) DEFAULT NULL,
  `up_date` varchar(500) DEFAULT NULL,
  `up_family` varchar(500) DEFAULT NULL,
  `up_taxo_lin` varchar(2000) DEFAULT NULL,
  `up_taxo_id` varchar(500) DEFAULT NULL,
  `up_cross_db_ref` varchar(2000) DEFAULT NULL,
  `binner` varchar(20) DEFAULT NULL,
  `min_contig_len` int(11) DEFAULT NULL,
  `k_mer_len` int(11) DEFAULT NULL,
  `marker` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=2 DEFAULT CHARSET=latin1;

CREATE TABLE `samples` (
  `id` int(11) NOT NULL,
  `samples_name` varchar(300) DEFAULT NULL,
  `samples_condition` varchar(300) DEFAULT NULL,
  `mg_sample1` varchar(300) DEFAULT NULL,
  `mg_sample2` varchar(300) DEFAULT NULL,
  `mg_description` text,
  `mt_sample1` varchar(300) DEFAULT NULL,
  `mt_sample2` varchar(300) DEFAULT NULL,
  `mt_description` text,
  `s_id` int(11) NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`s_id`)
) ENGINE=InnoDB AUTO_INCREMENT=16 DEFAULT CHARSET=latin1;

CREATE TABLE `exe_projects` (
  `exe_id` int(11) NOT NULL AUTO_INCREMENT,
  `samples` varchar(10000) DEFAULT NULL,
  `id` int(11) DEFAULT NULL,
  `s_id` int(11) DEFAULT NULL,
  PRIMARY KEY (`exe_id`)
) ENGINE=InnoDB AUTO_INCREMENT=2 DEFAULT CHARSET=latin1;

CREATE TABLE `users` (
 `id` int(11) NOT NULL AUTO_INCREMENT,
 `name` varchar(100) DEFAULT NULL,
 `username` varchar(30) DEFAULT NULL,
 `email` varchar(150) DEFAULT NULL,
 `password` varchar(100) DEFAULT NULL,
 `register_date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
 PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=5 DEFAULT CHARSET=latin1;


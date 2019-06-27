from flask import Flask, render_template, flash, redirect, request, url_for, session, logging, jsonify, send_file
#from data import Articles
from flask_mysqldb import MySQL
from wtforms import Form, widgets, BooleanField, StringField, TextAreaField, PasswordField, validators, FileField, SelectField, IntegerField, RadioField, SelectMultipleField
from passlib.hash import sha256_crypt
from functools import wraps
from werkzeug.utils import secure_filename
import os
import ast
import subprocess
from subprocess import Popen, PIPE, capture_output
from subprocess import check_output
import time
from ast import literal_eval
import glob

UPLOAD_FOLDER = '/home/mario/shared_folder/MOSCA_app/downloads'
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'])
app = Flask(__name__)
#images = Images(app)
app.secret_key='s$+adpowdadaopw31249ojfaés2i312948ecret123'
#app.secret_key='monkey'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
#app.root_path = ['/static']




#Config flask_mysqldb
app.config['MYSQL_HOST']='localhost'
app.config['MYSQL_USER']='root'
app.config['MYSQL_PASSWORD']=''
app.config['MYSQL_DB']='myproject'
app.config['MYSQL_CURSORCLASS']='DictCursor'
#init mySQL
mysql = MySQL(app)

state = 'True'
#Articles = Articles()

@app.route('/')
def index():
    return render_template('home.html')
#index
@app.route('/about')
def about():
    return render_template('about.html')

# Articles
@app.route('/articles')
def articles():
    # Create cursor
    cur = mysql.connection.cursor()

    # Get articles
    result = cur.execute("SELECT * FROM articles")

    articles = cur.fetchall()

    if result > 0:
        return render_template('articles.html', articles=articles)
    else:
        msg = 'No Projects Found'
        return render_template('articles.html', msg=msg)
    # Close connection
    cur.close()


#display articles
@app.route('/article/<string:id>/')
#def article(publish_id):
    #return render_template('article.html', publish_id = publish_id)
#display articles
def article(id):
    cur = mysql.connection.cursor()

    # Get article
    result = cur.execute("SELECT * FROM articles WHERE id = %s", [id])

    article = cur.fetchone()

    return render_template('article.html', article=article)


#register form class
class RegistrationForm(Form):
    name = StringField('Name', [validators.Length(min=4, max=50)])
    username = StringField('Username', [validators.Length(min=4, max=25)])
    email = StringField('Email Address', [validators.Length(min=6, max=35)])
    password = PasswordField('New Password', [
        validators.DataRequired(),
        validators.EqualTo('confirm', message='Passwords must match')
    ])
    confirm = PasswordField('Repeat Password')
    accept_tos = BooleanField('I accept the TOS', [validators.DataRequired()])

#user registration
@app.route('/register', methods=['GET', 'POST'])
def register():
    form = RegistrationForm(request.form)
    if request.method == 'POST' and form.validate():
        name = form.name.data
        email = form.username.data
        username = form.email.data
        password = sha256_crypt.encrypt(str(form.password.data))
        #user = User(form.name.data, form.username.data,
        #            form.email.data, form.password.data)
        #db_session.add(user)

        #CREATE CURSOR
        cur = mysql.connection.cursor()
        #execute
        cur.execute('INSERT INTO users(name, email, username, password) VALUES(%s, %s, %s, %s)',(name, email, username, password))

        # COMIT to db
        mysql.connection.commit()
        #close connection
        cur.close()
        flash('Thanks for registering you can now login', 'success')
        return redirect(url_for('index'))
    return render_template('register.html', form=form)

#user logging
@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        # Get Form Fields
        username = request.form['username']
        password_candidate = request.form['password']

        # Create cursor
        cur = mysql.connection.cursor()

        # Get user by username
        result = cur.execute("SELECT * FROM users WHERE username = %s", [username])

        if result > 0:
            # Get stored hash
            data = cur.fetchone()
            password = data['password']

            # Compare Passwords
            if sha256_crypt.verify(password_candidate, password):
                # Passed
                session['logged_in'] = True
                session['username'] = username

                flash('You are now logged in', 'success')
                return redirect(url_for('dashboard'))
            else:
                error = 'Invalid login'
                return render_template('login.html', error=error)
            # Close connection
            cur.close()
        else:
            error = 'Username not found'
            return render_template('login.html', error=error)

    return render_template('login.html')
#check if user is logged in:
def is_logged_in(f):
    @wraps(f)
    def wrap(*args, **kwargs):
        if 'logged_in' in session:
            return f(*args, **kwargs)
        else:
            flash('Unauthorized, Please login', 'danger')
            return redirect(url_for('login'))
    return wrap

#logout
@app.route('/logout')
@is_logged_in
def logout():
    session.clear()
    flash('You are now logged out', 'success')
    return redirect(url_for('login'))


#####################################################################################################################################################


#Create dashboard
@app.route('/dashboard')
@is_logged_in
def dashboard():
    global state
    state = 'True'
    #print(get_filelist())
    # Create cursor
    cur = mysql.connection.cursor()

    # Get articles
    #result = cur.execute("SELECT * FROM articles")
    # Show articles only from the user logged in
    result = cur.execute("SELECT * FROM projects")

    projects = cur.fetchall()
    #print(projects)
    projects_tags =  []
    for p in range(len(projects)):
        name = ''
        for names in projects[p]['author'].split(' '):
            name += names[0]
        project_tag = str(projects[p]['date']).split(' ')[0]+'_'+name+'_'+projects[p]['project_name']
        projects[p]['description'] = projects[p]['description'][3:-6]
        projects[p]['tags'] = project_tag

    if result > 0:
        return render_template('dashboard.html', projects=projects)
    else:
        msg = 'No Projects Found'
        cur.execute("ALTER TABLE projects AUTO_INCREMENT = 1")
        return render_template('dashboard.html', msg=msg)

    # Close connection
    cur.close()


#####################################################################################################################################################
#function to get files
def get_filelist(path):
    #path = '/media/sf_shared_folder/MOSCA_app/templates'
    files = []
    files.append(('','None'))
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            files.append((os.path.join(r, file),os.path.join(file)))
    return files

def get_filelist2(path):
    #path = '/media/sf_shared_folder/MOSCA_app/templates'
    files = []
    #files.append(('','None'))
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            files.append((os.path.join(r, file),os.path.join(file)))
    return files

#functions for directories
def get_dirlist(rootdir):

    dirlist = []

    with os.scandir(rootdir) as rit:
        for entry in rit:
            if not entry.name.startswith('.') and entry.is_dir():
                dirlist.append(entry.path)

    dirlist.sort() # Sort directory names
    return dirlist

#for multi check  box
def select_multi_checkbox(field, ul_class='', **kwargs):
    kwargs.setdefault('type', 'checkbox')
    field_id = kwargs.pop('id', field.id)
    html = [u'<ul %s>' % widgets.html_params(id=field_id, class_=ul_class)]
    for value, label, checked in field.iter_choices():
        choice_id = u'%s-%s' % (field_id, value)
        options = dict(kwargs, name=field.name, value=value, id=choice_id)
        if checked:
            options['checked'] = 'checked'

        if label in ['Entry name', 'Gene names', 'Protein names', 'EC number', 'Function[CC]', 'Pathway', 'Keywords', 'Protein existence', 'Gene ontology (GO)', 'Protein families', 'Taxonomic lineage (SUPERKINGDOM)',
            'Taxonomic lineage (PHYLUM)', 'Taxonomic lineage (CLASS)',
            'Taxonomic lineage (ORDER)', 'Taxonomic lineage (FAMILY)',
            'Taxonomic lineage (GENUS)', 'Taxonomic lineage (SPECIES)', 'BioCyc Collection of Pathway/Genome Databases',
             'BRENDA Comprehensive Enzyme Information System',
             'Conserved Domains Database',
             'evolutionary genealogy of genes: Non-supervised Orthologous Groups',
             'Ensembl eukaryotic genome annotation project',
             'Integrated resource of protein families, domains and functional sites',
             'KEGG: Kyoto Encyclopedia of Genes and Genomes',
             'KEGG Orthology (KO)', 'Pfam protein domain database',
             'Reactome - a knowledgebase of biological pathways and processes',
             'NCBI Reference Sequences',
             'UniPathway: a resource for the exploration and annotation of metabolic pathways']:
        #print(value,label,checked)
            checked = True
            html.append(u'<li><input %s checked/> ' % widgets.html_params(**options))
            html.append(u'<label for="%s">%s</label></li>' % (field_id, label))
        else:
            html.append(u'<li><input %s /> ' % widgets.html_params(**options))
            html.append(u'<label for="%s">%s</label></li>' % (field_id, label))
    html.append(u'</ul>')
    return u''.join(html)


class ProjectForm(Form):

    project_name = StringField('Project name', [validators.Length(min=1, max=300)])
    author = StringField('Author (1st and last name)', [validators.Length(min=1, max=300)])
    description = TextAreaField('Project description', [validators.Length(min=1,max=300)])

    #OUTPUT
    directories = get_dirlist('/media')
    drop_out_dir = []
    for i in directories:
        drop_out_dir.append((i,i))

    database_options = get_filelist2(os.path.join(app.instance_path)[0:-8]+'MOSCA/databases')
    #database_options = get_filelist2('/mnt/HDDSstorage/jsequeira/MOSCA/Databases/annotation_databases') ###path in server



    #database_dir = SelectField('Database file for Annotation', choices = database_options[::-1], default = database_options[0][1])
    database_dir = SelectField('Database file for Annotation ', default=[('/mosca/Databases/annotation_databases/uniprot.fasta','uniprot.fasta')],choices = [('/mosca/Databases/annotation_databases/uniprot.fasta','uniprot.fasta'),('/mosca/Databases/annotation_databases/ncbi.fasta','ncbi.fasta')])
    #database_dir = StringField('rRNA database directory', [validators.Length(min=1,max=300)])

    threads = IntegerField('Number of Threads', default=4)
    #threads = SelectField(choices = [('1','1 thread'),('2','2 threads'),('3','3 threads'),('4','4 threads')])
    sequencing = SelectField('Type of Sequencing', choices = [('PE','Paired-End'),('SE','Single-End')])
    quality_scores = SelectField(choices = [('phred33','phred33'),('phred64','phred64')])
    output_lvl = SelectField('Output Level', choices = [('minimum','minimum'),('medium','medium'),('maximum','maximum')])
    data_type = SelectField('Type of data coupled with metagenomics', choices = [('metatranscriptomics','metatranscriptomics'),('metaproteomics','metaproteomics')])

    preprocessing = BooleanField('Preprocessing',render_kw={'checked': True})
    assembly = BooleanField('Assembly',render_kw={'checked': True})
    binning = BooleanField('Binning',render_kw={'checked': True})

    #preprocessing = SelectField('No preprocessing?',choices = [(True,'Yes'),(False,'No')])
    #assembly = SelectField('No assembly?',choices = [(True,'Yes'),(False,'No')])
    #binning = SelectField('No binning?',choices = [(True,'Yes'),(False,'No')])

    #ASSEMBLY
    assembler = SelectField('Assembler', choices = [('metaSPAdes','metaSPAdes'),('MEGAHIT','MEGAHIT')])
    memory = IntegerField('Memory (in bytes) for assembly',validators=[validators.Optional()])
    k_mer_sizes = StringField('K mer sizes -  List of kmers for assembly (separated by comma), must be odd numbers. If empty default values will be used', validators=[validators.Optional()])
    #default for MegaHIT = [21,29,39,59,79,99,119,141]

    #MetaQUAST
    m = IntegerField('Minimum contig length to be considered for quality analysis of assembly', default = 500)

    #Bowtie2
    alignment_method = SelectField(choices = [('end-to-end','end-to-end'),('local','local')])
    alignment_options = SelectField(choices = [('sensitive','sensitive'),('very-fast','very-fast'),('fast','fast'),('very-sensitive','very-sensitive')])#adicionar local nas opções se o método for local


    #FragGeneScan
    #train_dir = SelectField('Train directory',choices = drop_out_dir)
    train = SelectField(choices = [('sanger_5','sanger_5'),('sanger_10','sanger_10'),('454_10','454_10'),('454_30','454_30'),('illumina_5','illumina_5'),('illumina_10','illumina_10')])

    #UniProt
    up_names_tax = SelectMultipleField('Names & Taxonomy', default=[('Entry name','Entry name')], choices = [('Entry name','Entry name'),('Gene names','Gene names'),('Gene names (ordered locus)','Gene names (ordered locus)'),('Gene names (ORF)','Gene names (ORF)'),('Gene names (primary)','Gene names (primary)'),('Gene names (synonym)','Gene names (synonym)'),('Organism','Organism'),('Organism ID','Organism ID'),('Protein names','Protein names'),('Proteomes','Proteomes'),('Taxonomic lineagea','Taxonomic lineage'),('Virus hosts','Virus hosts')], widget=select_multi_checkbox)

    up_sequences = SelectMultipleField('Sequences', choices = [('Alternative products (isoforms)','Alternative products (isoforms)'),('Alternative sequence','Alternative sequence'),('Erroneous gene model prediction','Erroneous gene model prediction'),('Fragment','Fragment'),('Gene encoded by','Gene encoded by'),('Length','Length'),('Mass','Mass'),('Mass spectrometry','Mass spectrometry'),('Natural variant','Natural variant'),('Non-adjacent residues','Non-adjacent residues'),('Non-standard residue','Non-standard residue'),('Non-terminal residue','Non-terminal residue'),('Polymorphism','Polymorphism'),('RNA editing','RNA editing'),('Sequence','Sequence'),('Sequence caution','Sequence caution'),('Sequence conflict','Sequence conflict'),('Sequence uncertainty','Sequence uncertainty'),('Sequence version','Sequence version')], widget=select_multi_checkbox)

    up_function = SelectMultipleField('Function', choices = [('Absorption','Absorption'),('Active site','Active site'),('Activity regulation','Activity regulation'),('Binding site','Binding site'),('Calcium binding','Calcium binding'),('Catalytic activity','Catalytic activity'),('Cofactor','Cofactor'),('DNA binding','DNA binding'),('EC number','EC number'),('Function[CC]','Function[CC]'),('Kinetics','Kinetics'),('Metal binding','Metal binding'),('Nucleotide binding','Nucleotide binding'),('Pathway','Pathway'),('pH dependence','pH dependence'),('Redox potential','Redox potential'),('Rhea Ids','Rhea Ids'),('Site','Site'),('Temperature dependence','Temperature dependence')], widget=select_multi_checkbox)

    up_miscellaneous =  SelectMultipleField('Miscellaneous', choices = [('Annotation','Annotation'),('Caution','Caution'),('Features','Features'),('Keyword ID','Keyword ID'),('Keywords','Keywords'),('Matched text','Matched text'),('Miscellaneous [CC]','Miscellaneous [CC]'),('Protein existence','Protein existence'),('Reviewed','Reviewed'),('Tools','Tools'),('UniParc','UniParc')],  widget=select_multi_checkbox)

    up_interaction = SelectMultipleField('Interaction', choices=[('Interacts with','Interacts with'),('Subunit structure [CC]','Subunit structure [CC]')],   widget=select_multi_checkbox)

    up_expression = SelectMultipleField('Expression',choices=[('Developmental stage','Developmental stage'),('Induction','Induction'),('Tissue specificity','Tissue specificity')], widget=select_multi_checkbox)

    up_gene_ont = SelectMultipleField('Gene Ontology (GO)',choices=[('Gene ontology (biological process)','Gene ontology (biological process)'),('Gene ontology (cellular component)','Gene ontology (cellular component)'),('Gene ontology (GO)','Gene ontology (GO)'),('Gene ontology (molecular function)','Gene ontology (molecular function)'),('Gene ontology IDs','Gene ontology IDs')], widget=select_multi_checkbox)

    up_chebi = SelectMultipleField('Chemical entities (ChEBI)', choices=[('ChEBI','ChEBI'),('ChEBI (Catalytic activity)','ChEBI (Catalytic activity)'),('ChEBI (Cofactor)','ChEBI (Cofactor)'),('ChEBI IDs','ChEBI IDs')], widget=select_multi_checkbox)

    up_path_biot = SelectMultipleField('Pathology & biotech', choices=[('Allergenic properties','Allergenic properties'),('Biotechnological use','Biotechnological use'),('Disruption phenotype','Disruption phenotype'),('Involvement in disease','Involvement in disease'),('Mutagenesis','Mutagenesis'),('Pharmaceutical use','Pharmaceutical use'),('Toxic dose','Toxic dose')], widget=select_multi_checkbox)

    up_cell_loc = SelectMultipleField('Subcellular location', choices=[('Intramembrane','Intramembrane'),('Subcellular location [CC]','Subcellular location [CC]'),('Topological domain','Topological domain'),('Transmembrane','Transmembrane')], widget=select_multi_checkbox)

    up_ptm = SelectMultipleField('PTM/Processing', choices=[('Chain','Chain'),('Cross-link','Cross-link'),('Disulfide bond','Disulfide bond'),('Glycosylation','Glycosylation'),('Initiator methionine','Initiator methionine'),('Lipidation','Lipidation'),('Modified residue','Modified residue'),('Peptide','Peptide'),('Post-translational modification','Post-translational modification'),('Propeptide','Propeptide'),('Signal peptide','Signal peptide'),('Transit peptide','Transit peptide')],  widget=select_multi_checkbox)

    up_structure = SelectMultipleField('Structure', choices=[('3D','3D'),('Beta strand','Beta strand'),('Helix','Helix'),('Turn','Turn')], widget=select_multi_checkbox)

    up_pubs = SelectMultipleField('Publications', choices=[('Mapped PubMed ID','Mapped PubMed ID'),('PubMed ID','PubMed ID')], widget=select_multi_checkbox)

    up_date = SelectMultipleField('Date of', choices=[('Date of creation','Date of creation'),('Date of last modification','Date of last modification'),('Date of last sequence modification','Date of last sequence modification'),('Entry version','Entry version')], widget=select_multi_checkbox)

    up_family = SelectMultipleField('Family & Domains', choices=[('Coiled coil','Coiled coil'),('Compositional bias','Compositional bias'),('Domain [CC]','Domain [CC]'),('Domain [FT]','Domain [FT]'),('Motif','Motif'),('Protein families','Protein families'),('Region','Region'),('Repeat','Repeat'),('Sequence similarities','Sequence similarities'),('Zinc finger','Zinc finger')], widget=select_multi_checkbox)

    up_taxo_lin= SelectMultipleField('Taxonomic lineage', choices=[('Taxonomic lineage (all)','Taxonomic lineage (all)'),('Taxonomic lineage (CLASS)','Taxonomic lineage (CLASS)'),('Taxonomic lineage (COHORT)','Taxonomic lineage (COHORT)'),('Taxonomic lineage (FAMILY)','Taxonomic lineage (FAMILY)'),('Taxonomic lineage (FORMA)','Taxonomic lineage (FORMA)'),('Taxonomic lineage (GENUS)','Taxonomic lineage (GENUS)'),('Taxonomic lineage (INFRACLASS)','Taxonomic lineage (INFRACLASS)'),('Taxonomic lineage (INFRAORDER)','Taxonomic lineage (INFRAORDER)'),('Taxonomic lineage (KINGDOM)','Taxonomic lineage (KINGDOM)'),('Taxonomic lineage (ORDER)','Taxonomic lineage (ORDER)'),('Taxonomic lineage (PARVORDER)','Taxonomic lineage (PARVORDER)'),('Taxonomic lineage (PHYLUM)','Taxonomic lineage (PHYLUM)'),('Taxonomic lineage (SPECIES)','Taxonomic lineage (SPECIES)'),('Taxonomic lineage (SPECIES_GROUP)','Taxonomic lineage (SPECIES_GROUP)'),('Taxonomic lineage (SPECIES_SUBGROUP)','Taxonomic lineage (SPECIES_SUBGROUP)'),('Taxonomic lineage (SUBCLASS)','Taxonomic lineage (SUBCLASS)'),('Taxonomic lineage (SUBCOHORT)','Taxonomic lineage (SUBCOHORT)'),('Taxonomic lineage (SUBFAMILY)','Taxonomic lineage (SUBFAMILY)'),('Taxonomic lineage (SUBGENUS)','Taxonomic lineage (SUBGENUS)'),('Taxonomic lineage (SUBKINGDOM)','Taxonomic lineage (SUBKINGDOM)'),('Taxonomic lineage (SUBORDER)','Taxonomic lineage (SUBORDER)'),('Taxonomic lineage (SUBPHYLUM)','Taxonomic lineage (SUBPHYLUM)'),('Taxonomic lineage (SUBSPECIES)','Taxonomic lineage (SUBSPECIES)'),('Taxonomic lineage (SUBTRIBE)','Taxonomic lineage (SUBTRIBE)'),('Taxonomic lineage (SUPERCLASS)','Taxonomic lineage (SUPERCLASS)'),('Taxonomic lineage (SUPERFAMILY)','Taxonomic lineage (SUPERFAMILY)'),('Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (SUPERKINGDOM)'),('Taxonomic lineage (SUPERORDER)','Taxonomic lineage (SUPERORDER)'),('Taxonomic lineage (SUPERPHYLUM)','Taxonomic lineage (SUPERPHYLUM)'),('Taxonomic lineage (TRIBE)','Taxonomic lineage (TRIBE)'),('Taxonomic lineage (VARIETAS)','Taxonomic lineage (VARIETAS)')], widget=select_multi_checkbox)

    up_taxo_id = SelectMultipleField('Taxonomic Identifier', choices =[('Taxonomic lineage IDs','Taxonomic lineage IDs')], widget=select_multi_checkbox)

    up_cross_db_ref = SelectMultipleField('', choices =[('Allergome; a platform for allergen knowledge', 'Allergome; a platform for allergen knowledge'), ('ArachnoServer: Spider toxin database', 'ArachnoServer: Spider toxin database'), ('Arabidopsis Information Portal', 'Arabidopsis Information Portal'), ('Bgee dataBase for Gene Expression Evolution', 'Bgee dataBase for Gene Expression Evolution'), ('BindingDB database of measured binding affinities', 'BindingDB database of measured binding affinities'), ('BioCyc Collection of Pathway/Genome Databases', 'BioCyc Collection of Pathway/Genome Databases'), ('The Biological General Repository for Interaction Datasets (BioGrid)', 'The Biological General Repository for Interaction Datasets (BioGrid)'), ('BioMuta curated single-nucleotide variation and disease association database', 'BioMuta curated single-nucleotide variation and disease association database'), ('BRENDA Comprehensive Enzyme Information System', 'BRENDA Comprehensive Enzyme Information System'), ('CarbonylDB database of protein carbonylation sites', 'CarbonylDB database of protein carbonylation sites'), ('Carbohydrate-Active enZymes', 'Carbohydrate-Active enZymes'), ('The Consensus CDS (CCDS) project', 'The Consensus CDS (CCDS) project'), ('Conserved Domains Database', 'Conserved Domains Database'), ('Candida Genome Database', 'Candida Genome Database'), ('ChEMBL database of bioactive drug-like small molecules', 'ChEMBL database of bioactive drug-like small molecules'), ('ChiTaRS: a database of human, mouse and fruit fly chimeric transcripts and RNA-sequencing data', 'ChiTaRS: a database of human, mouse and fruit fly chimeric transcripts and RNA-sequencing data'), ('CollecTF database of bacterial transcription factor binding sites', 'CollecTF database of bacterial transcription factor binding sites'), ('ComplexPortal: manually curated resource of macromolecular complexes', 'ComplexPortal: manually curated resource of macromolecular complexes'), ('2-DE database at Universidad Complutense de Madrid', '2-DE database at Universidad Complutense de Madrid'), ('ConoServer: Cone snail toxin database', 'ConoServer: Cone snail toxin database'), ('CORUM comprehensive resource of mammalian protein complexes', 'CORUM comprehensive resource of mammalian protein complexes'), ('Comparative Toxicogenomics Database', 'Comparative Toxicogenomics Database'), ('Database of single nucleotide polymorphism', 'Database of single nucleotide polymorphism'), ('DNA Data Bank of Japan; a nucleotide sequence database', 'DNA Data Bank of Japan; a nucleotide sequence database'), ('DEPOD human dephosphorylation database', 'DEPOD human dephosphorylation database'), ('Dictyostelium discoideum online informatics resource', 'Dictyostelium discoideum online informatics resource'), ('Database of interacting proteins', 'Database of interacting proteins'), ('DisGeNET', 'DisGeNET'), ('Database of protein disorder', 'Database of protein disorder'), ('Domain mapping of disease mutations (DMDM)', 'Domain mapping of disease mutations (DMDM)'), ('The DNASU plasmid repository', 'The DNASU plasmid repository'), ('DOSAC-COBS 2D-PAGE database', 'DOSAC-COBS 2D-PAGE database'), ('Drug and drug target database', 'Drug and drug target database'), ('EchoBASE - an integrated post-genomic database for E. coli', 'EchoBASE - an integrated post-genomic database for E. coli'), ('Escherichia coli strain K12 genome database', 'Escherichia coli strain K12 genome database'), ('evolutionary genealogy of genes: Non-supervised Orthologous Groups', 'evolutionary genealogy of genes: Non-supervised Orthologous Groups'), ('The Eukaryotic Linear Motif resource for Functional Sites in Proteins', 'The Eukaryotic Linear Motif resource for Functional Sites in Proteins'), ('EMBL nucleotide sequence database', 'EMBL nucleotide sequence database'), ('Ensembl eukaryotic genome annotation project', 'Ensembl eukaryotic genome annotation project'), ('Ensembl bacterial and archaeal genome annotation project', 'Ensembl bacterial and archaeal genome annotation project'), ('Ensembl fungal genome annotation project', 'Ensembl fungal genome annotation project'), ('Ensembl metazoan genome annotation project', 'Ensembl metazoan genome annotation project'), ('Ensembl plant genome annotation project', 'Ensembl plant genome annotation project'), ('Ensembl protists genome annotation project', 'Ensembl protists genome annotation project'), ('Enzyme nomenclature database', 'Enzyme nomenclature database'), ('Encyclopedia of Proteome Dynamics', 'Encyclopedia of Proteome Dynamics'), ('ESTHER database of the Alpha/Beta-hydrolase fold superfamily of proteins', 'ESTHER database of the Alpha/Beta-hydrolase fold superfamily of proteins'), ('European Hepatitis C Virus Database', 'European Hepatitis C Virus Database'), ('Eukaryotic Pathogen Database Resources', 'Eukaryotic Pathogen Database Resources'), ('Relative evolutionary importance of amino acids within a protein sequence', 'Relative evolutionary importance of amino acids within a protein sequence'), ('ExpressionAtlas, Differential and Baseline Expression', 'ExpressionAtlas, Differential and Baseline Expression'), ('Drosophila genome database', 'Drosophila genome database'), ('GenAtlas: human gene database', 'GenAtlas: human gene database'), ('GenBank nucleotide sequence database', 'GenBank nucleotide sequence database'), ('Gene3D Structural and Functional Annotation of Protein Families', 'Gene3D Structural and Functional Annotation of Protein Families'), ('GeneCards: human genes, protein and diseases', 'GeneCards: human genes, protein and diseases'), ('GeneDB pathogen genome database from Sanger Institute', 'GeneDB pathogen genome database from Sanger Institute'), ('Database of genes from NCBI RefSeq genomes', 'Database of genes from NCBI RefSeq genomes'), ('GeneReviews a resource of expert-authored, peer-reviewed disease descriptions.', 'GeneReviews a resource of expert-authored, peer-reviewed disease descriptions.'), ('Ensembl GeneTree', 'Ensembl GeneTree'), ('Genevisible search portal to normalized and curated expression data from Genevestigator', 'Genevisible search portal to normalized and curated expression data from Genevestigator'), ('The Gene Wiki collection of pages on human genes and proteins', 'The Gene Wiki collection of pages on human genes and proteins'), ('Database of phenotypes from RNA interference screens in Drosophila and Homo sapiens', 'Database of phenotypes from RNA interference screens in Drosophila and Homo sapiens'), ('GlyConnect protein glycosylation platform', 'GlyConnect protein glycosylation platform'), ('Gene Ontology', 'Gene Ontology'), ('Information system for G protein-coupled receptors (GPCRs)', 'Information system for G protein-coupled receptors (GPCRs)'), ('Gramene; a comparative resource for plants', 'Gramene; a comparative resource for plants'), ('IUPHAR/BPS Guide to PHARMACOLOGY', 'IUPHAR/BPS Guide to PHARMACOLOGY'), ('H-Invitational Database, human transcriptome db', 'H-Invitational Database, human transcriptome db'), ('HAMAP database of protein families', 'HAMAP database of protein families'), ('Human Gene Nomenclature Database', 'Human Gene Nomenclature Database'), ('The HOGENOM Database of Homologous Genes from Fully Sequenced Organisms', 'The HOGENOM Database of Homologous Genes from Fully Sequenced Organisms'), ('Human Protein Atlas', 'Human Protein Atlas'), ('Human Unidentified Gene-Encoded large proteins database', 'Human Unidentified Gene-Encoded large proteins database'), ('The international ImMunoGeneTics information system', 'The international ImMunoGeneTics information system'), ('InParanoid: Eukaryotic Ortholog Groups', 'InParanoid: Eukaryotic Ortholog Groups'), ('Protein interaction database and analysis system', 'Protein interaction database and analysis system'), ('Integrated resource of protein families, domains and functional sites', 'Integrated resource of protein families, domains and functional sites'), ('iPTMnet integrated resource for PTMs in systems biology context', 'iPTMnet integrated resource for PTMs in systems biology context'), ('jPOST - Japan Proteome Standard Repository/Database', 'jPOST - Japan Proteome Standard Repository/Database'), ('KEGG: Kyoto Encyclopedia of Genes and Genomes', 'KEGG: Kyoto Encyclopedia of Genes and Genomes'), ('KEGG Orthology (KO)', 'KEGG Orthology (KO)'), ('Legionella pneumophila genome database', 'Legionella pneumophila genome database'), ('Mycobacterium leprae genome database', 'Mycobacterium leprae genome database'), ('Maize Genetics and Genomics Database', 'Maize Genetics and Genomics Database'), ('MalaCards human disease database', 'MalaCards human disease database'), ('MaxQB - The MaxQuant DataBase', 'MaxQB - The MaxQuant DataBase'), ('MEROPS protease database', 'MEROPS protease database'), ('Mouse genome database (MGD) from Mouse Genome Informatics (MGI)', 'Mouse genome database (MGD) from Mouse Genome Informatics (MGI)'), ('Microbial advanced database', 'Microbial advanced database'), ('Online Mendelian Inheritance in Man (OMIM)', 'Online Mendelian Inheritance in Man (OMIM)'), ('Molecular INTeraction database', 'Molecular INTeraction database'), ('MobiDB: a database of protein disorder and mobility annotations', 'MobiDB: a database of protein disorder and mobility annotations'), ('Database of comparative protein structure models', 'Database of comparative protein structure models'), ('MoonDB Database of extreme multifunctional and moonlighting proteins', 'MoonDB Database of extreme multifunctional and moonlighting proteins'), ('MoonProt database of moonlighting proteins', 'MoonProt database of moonlighting proteins'), ('mycoCLAP, a database of fungal genes encoding lignocellulose-active proteins', 'mycoCLAP, a database of fungal genes encoding lignocellulose-active proteins'), ('neXtProt; the human protein knowledge platform', 'neXtProt; the human protein knowledge platform'), ('USC-OGP 2-DE database', 'USC-OGP 2-DE database'), ('Identification of Orthologs from Complete Genome Data', 'Identification of Orthologs from Complete Genome Data'), ('Open Targets', 'Open Targets'), ('Orphanet; a database dedicated to information on rare diseases and orphan drugs', 'Orphanet; a database dedicated to information on rare diseases and orphan drugs'), ('Database of Orthologous Groups', 'Database of Orthologous Groups'), ('The PANTHER Classification System', 'The PANTHER Classification System'), ('Pathosystems Resource Integration Center (PATRIC)', 'Pathosystems Resource Integration Center (PATRIC)'), ('PaxDb, a database of protein abundance averages across all three domains of life', 'PaxDb, a database of protein abundance averages across all three domains of life'), ('Protein Data Bank Europe', 'Protein Data Bank Europe'), ('Protein Data Bank Japan', 'Protein Data Bank Japan'), ('PDBsum; at-a-glance overview of macromolecular structures', 'PDBsum; at-a-glance overview of macromolecular structures'), ('PeptideAtlas', 'PeptideAtlas'), ('PeroxiBase, a peroxidase database', 'PeroxiBase, a peroxidase database'), ('Pfam protein domain database', 'Pfam protein domain database'), ('The Pharmacogenetics and Pharmacogenomics Knowledge Base', 'The Pharmacogenetics and Pharmacogenomics Knowledge Base'), ('Comprehensive resource for the study of protein post-translational modifications (PTMs) in human, mouse and rat.', 'Comprehensive resource for the study of protein post-translational modifications (PTMs) in human, mouse and rat.'), ('Database for complete collections of gene phylogenies', 'Database for complete collections of gene phylogenies'), ('Protein sequence database of the Protein Information Resource', 'Protein sequence database of the Protein Information Resource'), ('PIRSF; a whole-protein classification database', 'PIRSF; a whole-protein classification database'), ('CutDB - Proteolytic event database', 'CutDB - Proteolytic event database'), ('Schizosaccharomyces pombe database', 'Schizosaccharomyces pombe database'), ('PRoteomics IDEntifications database', 'PRoteomics IDEntifications database'), ('Protein Motif fingerprint database; a protein domain database', 'Protein Motif fingerprint database; a protein domain database'), ('Protein Ontology', 'Protein Ontology'), ('ProDom; a protein domain database', 'ProDom; a protein domain database'), ('Protein Mass spectra EXtraction', 'Protein Mass spectra EXtraction'), ('PROSITE; a protein domain and family database', 'PROSITE; a protein domain and family database'), ('Proteomes', 'Proteomes'), ('ProteomicsDB human proteome resource', 'ProteomicsDB human proteome resource'), ('ProtoNet; Automatic hierarchical classification of proteins', 'ProtoNet; Automatic hierarchical classification of proteins'), ('Pseudomonas genome database', 'Pseudomonas genome database'), ('Protein Data Bank RCSB', 'Protein Data Bank RCSB'), ('Reactome - a knowledgebase of biological pathways and processes', 'Reactome - a knowledgebase of biological pathways and processes'), ('Restriction enzymes and methylases database', 'Restriction enzymes and methylases database'), ('NCBI Reference Sequences', 'NCBI Reference Sequences'), ('REPRODUCTION-2DPAGE', 'REPRODUCTION-2DPAGE'), ('Rat genome database', 'Rat genome database'), ('Rodent Unidentified Gene-Encoded large proteins database', 'Rodent Unidentified Gene-Encoded large proteins database'), ('SABIO-RK: Biochemical Reaction Kinetics Database', 'SABIO-RK: Biochemical Reaction Kinetics Database'), ('The Structural Biology Knowledgebase', 'The Structural Biology Knowledgebase'), ('Structure-Function Linkage Database', 'Structure-Function Linkage Database'), ('Saccharomyces Genome Database', 'Saccharomyces Genome Database'), ('SignaLink: a signaling pathway resource with multi-layered regulatory networks', 'SignaLink: a signaling pathway resource with multi-layered regulatory networks'), ('SIGNOR Signaling Network Open Resource', 'SIGNOR Signaling Network Open Resource'), ('Simple Modular Architecture Research Tool; a protein domain database', 'Simple Modular Architecture Research Tool; a protein domain database'), ('SWISS-MODEL Repository - a database of annotated 3D protein structure models', 'SWISS-MODEL Repository - a database of annotated 3D protein structure models'), ('The Stanford Online Universal Resource for Clones and ESTs', 'The Stanford Online Universal Resource for Clones and ESTs'), ('STRING: functional protein association networks', 'STRING: functional protein association networks'), ('Superfamily database of structural and functional annotation', 'Superfamily database of structural and functional annotation'), ('Two-dimensional polyacrylamide gel electrophoresis database from the Geneva University Hospital', 'Two-dimensional polyacrylamide gel electrophoresis database from the Geneva University Hospital'), ('SWISS-MODEL Interactive Workspace', 'SWISS-MODEL Interactive Workspace'), ('SwissLipids knowledge resource for lipid biology', 'SwissLipids knowledge resource for lipid biology'), ('SwissPalm database of S-palmitoylation events', 'SwissPalm database of S-palmitoylation events'), ('The Arabidopsis Information Resource', 'The Arabidopsis Information Resource'), ('Transport Classification Database', 'Transport Classification Database'), ('TIGRFAMs; a protein family database', 'TIGRFAMs; a protein family database'), ('Consortium for Top Down Proteomics', 'Consortium for Top Down Proteomics'), ('TreeFam database of animal gene trees', 'TreeFam database of animal gene trees'), ('Mycobacterium tuberculosis strain H37Rv genome database', 'Mycobacterium tuberculosis strain H37Rv genome database'), ('University College Dublin 2-DE Proteome Database', 'University College Dublin 2-DE Proteome Database'), ('UCSC genome browser', 'UCSC genome browser'), ('UniCarbKB; an annotated and curated database of glycan structures', 'UniCarbKB; an annotated and curated database of glycan structures'), ('UniLectin database of carbohydrate-binding proteins', 'UniLectin database of carbohydrate-binding proteins'), ('UniPathway: a resource for the exploration and annotation of metabolic pathways', 'UniPathway: a resource for the exploration and annotation of metabolic pathways'), ('Bioinformatics Resource for Invertebrate Vectors of Human Pathogens', 'Bioinformatics Resource for Invertebrate Vectors of Human Pathogens'), ('Vertebrate Gene Nomenclature Database', 'Vertebrate Gene Nomenclature Database'), ('WormBase ParaSite', 'WormBase ParaSite'), ('The World-2DPAGE database', 'The World-2DPAGE database'), ('WormBase', 'WormBase'), ('Xenopus laevis and tropicalis biology and genomics resource', 'Xenopus laevis and tropicalis biology and genomics resource'), ('Zebrafish Information Network genome database', 'Zebrafish Information Network genome database')], widget=select_multi_checkbox)

    #VizBin
    binner = SelectField('Binner tool',default=[('VizBin','VizBin')],choices = [('VizBin','VizBin'),('MaxBin2','MaxBin2')])
    min_contig_len = IntegerField('Minimum contig length', default=1000)
    k_mer_len = IntegerField('k-mer length', default=5)
    marker = SelectField('Gene Marker set',choices = [('40 ','40'),('107','107')])

@app.route('/add_article', methods = ['GET', 'POST'])
@is_logged_in
def add_article():
    form = ProjectForm(request.form)
    if request.method == 'POST' and form.validate():
        project_name = form.project_name.data
        author = form.author.data
        description = form.description.data

        #OUTPUT
        database_dir = form.database_dir.data
        threads = form.threads.data
        sequencing = form.sequencing.data
        quality_scores = form.quality_scores.data
        output_lvl = form.output_lvl.data
        data_type = form.data_type.data

        preprocessing = form.preprocessing.data
        assembly = form.assembly.data
        binning = form.binning.data

        #ASSEMBLY
        assembler = form.assembler.data
        memory = form.memory.data
        k_mer_sizes = form.k_mer_sizes.data

        #MetaQUAST
        m = form.m.data

        #Bowtie2

        alignment_method = form.alignment_method.data
        alignment_options = form.alignment_options.data


        #FragGeneScan
        #train_dir = form.train_dir.data
        train = form.train.data

        #UniProt
        up_names_tax = form.up_names_tax.data
        up_sequences = form.up_sequences.data
        up_function = form.up_function.data
        up_miscellaneous = form.up_miscellaneous.data
        up_interaction = form.up_interaction.data
        up_expression = form.up_expression.data
        up_gene_ont = form.up_gene_ont.data
        up_chebi = form.up_chebi.data
        up_path_biot = form.up_path_biot.data
        up_cell_loc = form.up_cell_loc.data
        up_ptm = form.up_ptm.data
        up_structure = form.up_structure.data
        up_pubs = form.up_pubs.data
        up_date = form.up_date.data
        up_family = form.up_family.data
        up_taxo_lin = form.up_taxo_lin.data
        up_taxo_id = form.up_taxo_id.data
        up_cross_db_ref = form.up_cross_db_ref.data

        #VizBin
        binner = form.binner.data
        min_contig_len = form.min_contig_len.data
        k_mer_len = form.k_mer_len.data
        marker = form.marker.data

        #create CURSOR
        cur = mysql.connection.cursor()
        #execute
        if alignment_method == 'local':
            alignment_options = alignment_options + '-local'

        cur.execute('INSERT INTO projects(project_name,author,description,database_dir,threads,sequencing,quality_scores, output_lvl, data_type, preprocessing, assembly, binning,assembler,memory,k_mer_sizes,m,alignment_method,alignment_options, train, up_names_tax, up_sequences, up_function, up_miscellaneous,up_interaction,up_expression,up_gene_ont,up_chebi,up_path_biot,up_cell_loc,up_ptm,up_structure,up_pubs,up_date,up_family,up_taxo_lin,up_taxo_id,up_cross_db_ref,binner, min_contig_len, k_mer_len,marker) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (project_name,author,description,database_dir,threads,sequencing,quality_scores,output_lvl,data_type,preprocessing, assembly, binning,assembler,memory,k_mer_sizes,m,alignment_method,alignment_options, train, str(up_names_tax), str(up_sequences), str(up_function), str(up_miscellaneous),str(up_interaction),str(up_expression),str(up_gene_ont),str(up_chebi),str(up_path_biot),str(up_cell_loc),str(up_ptm),str(up_structure),str(up_pubs),str(up_date),str(up_family),str(up_taxo_lin),str(up_taxo_id),str(up_cross_db_ref),binner, min_contig_len, k_mer_len, marker))

        #commit to db
        mysql.connection.commit()
        #close connection
        cur.close()
        flash('Parameters added','success')
        return redirect(url_for('dashboard'))
    return render_template('add_article.html', form=form)

#####################################################################################################################################################


#SAMPLES

class SamplesForm(Form):
    #directories = get_filelist()
    drop_out_dir = get_filelist(os.path.join(app.instance_path)[0:-8]+'input_files')
    #print(drop_out_dir)
    #drop_out_dir.append(('','None'))
    #for i in directories:
    #    drop_out_dir.append((i,i))
    samples_name = StringField('Name your samples set:')
    samples_condition = StringField('Samples condition')
    #mg_path = SelectField('Metagenomic files directory' , choices = drop_out_dir)
    mg_sample1 = SelectField('Metagenomic file R1' , choices = drop_out_dir)
    mg_sample2 = SelectField('Metagenomic file R2', choices = drop_out_dir)
    mg_description = StringField('Metagenomic samples Description')

    #mt_path = SelectField('Metatranscriptomics files directory', choices = drop_out_dir)
    mt_sample1 = SelectField('Metatranscriptomic file R1', choices = drop_out_dir)
    mt_sample2 = SelectField('Metatranscriptomic file R2', choices = drop_out_dir)
    mt_description = StringField('Metatranscriptomic samples Description')


@app.route('/add_data/<string:id>', methods = ['GET', 'POST'])
@is_logged_in
def add_data(id):
    form = SamplesForm(request.form)
    cur = mysql.connection.cursor()
    result = cur.execute("SELECT * FROM samples WHERE id=%s",[id])
    samples = cur.fetchall()
    result2 = cur.execute("SELECT project_name,author,date FROM projects WHERE id=%s",[id])
    name = cur.fetchone()
    n = name['author'].split(' ')
    tag_n= ''
    for i in n:
        tag_n += i[0]
    name = name['project_name']+'_'+tag_n+'_'+str(name['date']).split(' ')[0]
    result3 = cur.execute("SELECT count(s_id) FROM samples WHERE id=%s",[id])
    n_samples = cur.fetchone()
    n_samples = n_samples['count(s_id)']
    if request.method == 'POST' and form.validate():
        samples_name = form.samples_name.data
        samples_condition = form.samples_condition.data
        #mg_path = form.mg_path.data
        mg_sample1 = form.mg_sample1.data
        mg_sample2 = form.mg_sample2.data
        mg_description = form.mg_description.data
        #mt_path = form.mt_path.data
        mt_sample1 = form.mt_sample1.data
        mt_sample2 = form.mt_sample2.data
        mt_description = form.mt_description.data

        #execute
        if mg_sample1 and not mt_sample1:
            if mg_sample2 and not mt_sample2:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_sample1,mg_sample2,mg_description,mt_description) VALUES(%s,%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition,mg_sample1,mg_sample2,
                mg_description, mt_description))
            elif not mg_sample2 and not mt_sample2:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_sample1,mg_description,mt_description) VALUES(%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition, mg_sample1,
                mg_description, mt_description))
        elif mt_sample1 and not mg_sample1:
            if mt_sample2 and not mg_sample2:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_description,mt_sample1,mt_sample2,mt_description) VALUES(%s,%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition,
                mg_description, mt_sample1, mt_sample2, mt_description))
            elif not mt_sample2 and not mg_sample2:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_description,mt_sample1,mt_description) VALUES(%s,%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition,
                mg_description, mt_sample1, mt_description))
        elif mg_sample1 and mg_sample1:
            if mg_sample2 and not mt_sample2:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_sample1,mg_sample2,mg_description,mt_sample1,mt_description) VALUES(%s,%s,%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition, mg_sample1, mg_sample2,
                mg_description, mt_sample1, mt_description))
            elif mt_sample2 and not mg_sample2:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_sample1,mg_description,mt_sample1,mt_sample2,mt_description) VALUES(%s,%s,%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition, mg_sample1,
                mg_description, mt_sample1, mt_sample2, mt_description))
            elif not mg_sample2 and not mt_sample2:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_sample1,mg_description,mt_sample1,mt_description) VALUES(%s,%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition, mg_sample1,
                mg_description, mt_sample1, mt_description))
            else:
                cur.execute('INSERT INTO samples(id,samples_name,samples_condition,mg_sample1,mg_sample2,mg_description,mt_sample1,mt_sample2,mt_description) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s)',(id,samples_name,samples_condition, mg_sample1, mg_sample2, mg_description, mt_sample1, mt_sample2, mt_description))
        else:
            flash('Correct samples input','danger')
            return redirect(request.url)




        #commit to db
        mysql.connection.commit()
        #close connection
        cur.close()
        flash('Parameters added','success')
        return redirect(url_for('add_data', id = id, name = name, n_samples=n_samples))

    if result > 0:
        return render_template('add_data.html', samples = samples, form=form, id = id, name = name, n_samples=n_samples)
    else:
        msg = 'No Samples Found'
        return render_template('add_data.html', msg=msg, form=form, id = id, name = name, n_samples=n_samples)

#####################################################################################################################################################


#EDIT DATA
@app.route('/edit_sample/<string:id>', methods = ['GET', 'POST'])
@is_logged_in
def edit_sample(id):

    form = SamplesForm(request.form)
    cur = mysql.connection.cursor()

    result = cur.execute('SELECT * FROM samples WHERE s_id = %s',[id])

    sample = cur.fetchone()
    s_id = sample['s_id']
    p_id = sample['id']
    s_name = sample['samples_name']
    file=sample['mg_sample1'].split('/')[-1]
    print(file)
    form.samples_name.data = sample['samples_name']
    form.samples_condition.data = sample['samples_condition']
    #form.mg_path.data = sample['mg_path']
    form.mg_sample1.data = sample['mg_sample1']
    form.mg_sample2.data = sample['mg_sample2']
    form.mg_description.data = sample['mg_description']
    #form.mt_path.data= sample['mt_path']
    form.mt_sample1.data = sample['mt_sample1']
    form.mt_sample2.data = sample['mt_sample2']
    form.mt_description.data = sample['mt_description']

    cur.close()


    if request.method == 'POST' :

        samples_name = request.form.get('samples_name')
        samples_condition = request.form.get('samples_condition')
        #mg_path = request.form.get('mg_path')
        mg_sample1 = request.form.get('mg_sample1')
        mg_sample2 = request.form.get('mg_sample2')
        mg_description = request.form.get('mg_description')
        #mt_path = request.form.get('mt_path')
        mt_sample1 = request.form.get('mt_sample1')
        mt_sample2 = request.form.get('mt_sample2')
        mt_description = request.form.get('mt_description')


        #create CURSOR
        cur = mysql.connection.cursor()
        #execute
        if mg_sample1 and not mt_sample1:
            if mg_sample2 and not mt_sample2:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition, mg_sample1, mg_sample2, mg_description,None, None, None, mt_description, s_id))
            elif not mg_sample2 and not mt_sample2:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition, mg_sample1, None, mg_description,None, None, None, mt_description, s_id))

        elif mt_sample1 and not mg_sample1:
            if mt_sample2 and not mg_sample2:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition,None, None, None, mg_description, mt_sample1, mt_sample2, mt_description, s_id))
            elif not mt_sample2 and not mg_sample2:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition,None, None, None, mg_description, mt_sample1, None, mt_description, s_id))

        elif mg_sample1 and mg_sample1:
            if mg_sample2 and not mt_sample2:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition, mg_sample1, mg_sample2, mg_description, mt_sample1, None, mt_description, s_id))
            elif mt_sample2 and not mg_sample2:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition, mg_sample1, None, mg_description, mt_sample1, mt_sample2, mt_description, s_id))
            elif not mg_sample2 and not mt_sample2:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition, mg_sample1, None, mg_description, mt_sample1, None, mt_description, s_id))
            else:
                cur.execute('UPDATE samples SET samples_name=%s, samples_condition=%s, mg_sample1=%s, mg_sample2=%s, mg_description=%s, mt_sample1=%s, mt_sample2=%s, mt_description=%s WHERE s_id=%s',(samples_name,samples_condition, mg_sample1, mg_sample2, mg_description, mt_sample1, mt_sample2, mt_description, s_id))

        #commit to db
        mysql.connection.commit()
        #close connection
        cur.close()
        flash('Sample Updated','success')
        return redirect(url_for('add_data', id = p_id))

    mysql.connection.commit()
    #close connection
    cur.close()
    return render_template('edit_sample.html', form=form, id = p_id, s_name=s_name, s_id = s_id)

#####################################################################################################################################################


#EDIT ARTICLE
#Precisa de ser atualizado
@app.route('/edit_article/<string:id>', methods = ['GET', 'POST'])
@is_logged_in
def edit_article(id):
    #create CURSOR
    cur = mysql.connection.cursor()

    result = cur.execute('SELECT * FROM projects WHERE id = %s',[id])

    project = cur.fetchone()


    cur.close()


    #get form

    form = ProjectForm(request.form)

    #populate article form Fields
    form.project_name.data = project['project_name']
    form.author.data = project['author']
    form.description.data = project['description']
    form.database_dir.data = project['database_dir']
    form.threads.data = project['threads']
    form.sequencing.data = project['sequencing']
    form.quality_scores.data = project['quality_scores']

    form.output_lvl.data = project['output_lvl']
    form.data_type.data = project['data_type']

    form.preprocessing.data = project['preprocessing']
    form.assembly.data = project['assembly']
    form.binning.data = project['binning']

    form.assembler.data = project['assembler']
    form.memory.data = project['memory']
    form.k_mer_sizes.data = project['k_mer_sizes']
    form.m.data = project['m']
    form.alignment_method.data = project['alignment_method']
    form.alignment_options.data = project['alignment_options']
    #form.train_dir.data = project['train_dir']
    form.train.data = project['train']
    form.up_names_tax.data = project['up_names_tax']
    form.up_sequences.data = project['up_sequences']
    form.up_function.data = project['up_function']
    form.up_miscellaneous.data = project['up_miscellaneous']
    form.up_interaction.data = project['up_interaction']
    form.up_expression.data = project['up_expression']
    form.up_gene_ont.data = project['up_gene_ont']
    form.up_chebi.data = project['up_chebi']
    form.up_path_biot.data = project['up_path_biot']
    form.up_cell_loc.data = project['up_cell_loc']
    form.up_ptm.data = project['up_ptm']
    form.up_structure.data = project['up_structure']
    form.up_pubs.data = project['up_pubs']
    form.up_date.data = project['up_date']
    form.up_family.data = project['up_family']
    form.up_taxo_lin.data = project['up_taxo_lin']
    form.up_taxo_id.data = project['up_taxo_id']
    form.up_cross_db_ref.data = project['up_cross_db_ref']
    form.binner.data = project['binner']
    form.min_contig_len.data = project['min_contig_len']
    form.k_mer_len.data = project['k_mer_len']
    form.marker.data = project['marker']

    print(form.up_names_tax.data)
    if request.method == 'POST' and form.validate():

        project_name = request.form['project_name']
        author = request.form['author']
        description = request.form['description']

        #OUTPUT
        database_dir = request.form['database_dir']
        threads = request.form['threads']
        sequencing = request.form['sequencing']
        quality_scores = request.form['quality_scores']
        output_lvl = request.form['output_lvl']
        data_type = request.form['data_type']
        #ASSEMBLY
        assembler = request.form['assembler']
        memory = request.form['memory']
        k_mer_sizes = request.form['k_mer_sizes']

        #MetaQUAST
        m = request.form['m']

        #Bowtie2

        alignment_method = request.form['alignment_method']
        alignment_options = request.form['alignment_options']


        #FragGeneScan
        #train_dir = request.form['train_dir']
        train = request.form['train']

        #UniProt
        up_names_tax = request.form['up_names_tax']
        up_sequences = request.form.get['up_sequences']
        up_function = request.form.get['up_function']
        up_miscellaneous = request.form['up_miscellaneous']
        up_interaction = request.form['up_interaction']
        up_expression = request.form['up_expression']
        up_gene_ont = request.form['up_gene_ont']
        up_chebi = request.form['up_chebi']
        up_path_biot = request.form['up_path_biot']
        up_cell_loc = request.form['up_cell_loc']
        up_ptm = request.form['up_ptm']
        up_structure = request.form['up_structure']
        up_pubs = request.form['up_pubs']
        up_date = request.form['up_date']
        up_family = request.form['up_family']
        up_taxo_lin = request.form['up_taxo_lin']
        up_taxo_id = request.form['up_taxo_id']
        up_cross_db_ref = request.form['up_cross_db_ref']

        #VizBin
        binner = request.form['binner']
        min_contig_len = request.form['min_contig_len']
        k_mer_len = request.form['k_mer_len']
        marker = request.form['marker']

        #create CURSOR
        cur = mysql.connection.cursor()
        #execute
        if alignment_method == 'local':
            alignment_options = alignment_options + '-local'

        #cur.execute('DELETE FROM projects WHERE id= %s',[id])

        #cur.execute('INSERT INTO projects(project_name,author,description,database_dir,threads,sequencing,quality_scores, output_lvl, data_type, preprocessing, assembly, binning,assembler,memory,k_mer_sizes,m,alignment_method,alignment_options, train, up_names_tax, up_sequences, up_function, up_miscellaneous,up_interaction,up_expression,up_gene_ont,up_chebi,up_path_biot,up_cell_loc,up_ptm,up_structure,up_pubs,up_date,up_family,up_taxo_lin,up_taxo_id,up_cross_db_ref,binner, min_contig_len, k_mer_len,marker) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', (project_name,author,description,database_dir,threads,sequencing,quality_scores,output_lvl,data_type,preprocessing, assembly, binning,assembler,memory,k_mer_sizes,m,alignment_method,alignment_options, train, str(up_names_tax), str(up_sequences), str(up_function), str(up_miscellaneous),str(up_interaction),str(up_expression),str(up_gene_ont),str(up_chebi),str(up_path_biot),str(up_cell_loc),str(up_ptm),str(up_structure),str(up_pubs),str(up_date),str(up_family),str(up_taxo_lin),str(up_taxo_id),str(up_cross_db_ref),binner, min_contig_len, k_mer_len, marker))


        cur.execute('UPDATE projects SET project_name=%s, author=%s, description=%s, database_dir=%s, threads=%s, sequencing=%s, quality_scores=%s, output_lvl=%s, data_type=%s, assembler=%s, memory=%s, k_mer_sizes=%s, m=%s, alignment_method=%s, alignment_options=%s, train=%s, up_names_tax=%s, up_sequences=%s, up_function=%s, up_miscellaneous=%s, up_interaction=%s, up_expression=%s, up_gene_ont=%s, up_chebi=%s, up_path_biot=%s, up_cell_loc=%s, up_ptm=%s, up_structure=%s, up_pubs=%s, up_date=%s, up_family=%s, up_taxo_lin=%s, up_taxo_id=%s, up_cross_db_ref=%s, binner=%s, min_contig_len=%s, k_mer_len=%s, marker=%s WHERE id=%s',(	project_name,author,description,database_dir,threads,sequencing,quality_scores,output_lvl,data_type,assembler,memory,k_mer_sizes,m,alignment_method,alignment_options,train,up_names_tax[0],up_sequences,up_function,up_miscellaneous,up_interaction,up_expression,up_gene_ont,up_chebi,up_path_biot,up_cell_loc,up_ptm,up_structure,up_pubs,up_date,up_family,up_taxo_lin,up_taxo_id,up_cross_db_ref,binner,min_contig_len,k_mer_len,marker))
        #commit to db
        mysql.connection.commit()
        #close connection
        cur.close()
        flash('Project Updated','success')
        return redirect(url_for('dashboard'))
    return render_template('edit_article.html', form=form)

#####################################################################################################################################################


#DELETE ARTICLE
@app.route('/delete_article/<string:id>', methods = ['POST'])
@is_logged_in
def delete_article(id):
    #create CURSOR
    cur = mysql.connection.cursor()
    #execute
    cur.execute('DELETE FROM projects WHERE id= %s',[id])
    cur.execute('DELETE FROM samples WHERE id= %s',[id])
    #commit to db
    mysql.connection.commit()
    #close connection
    cur.close()
    flash('Project Deleted','success')
    return redirect(url_for('dashboard'))

#####################################################################################################################################################

#DELETE SAMPLES FROM A PROJECT
@app.route('/delete_samples/<string:id1>/<string:id2>', methods = ['POST'])
@is_logged_in
def delete_samples(id1,id2):
    #create CURSOR
    cur = mysql.connection.cursor()
    #execute
    cur.execute('DELETE FROM samples WHERE s_id= %s',[id2])
    #commit to db
    mysql.connection.commit()
    #close connection
    cur.close()
    flash('Project Samples','success')
    return redirect(url_for('add_data', id= id1))


#####################################################################################################################################################


#for multi check  box
def select_multi_checkbox2(field, ul_class='', **kwargs):
    kwargs.setdefault('type', 'checkbox')
    field_id = kwargs.pop('id', field.id)
    html = [u'<ul %s>' % widgets.html_params(id=field_id, class_=ul_class)]
    for value, label, checked in field.iter_choices():
        choice_id = u'%s-%s' % (field_id, value)
        options = dict(kwargs, name=field.name, value=value, id=choice_id)
        if checked:
            options['checked'] = 'checked'

        if label == 'Entry name':
        #print(value,label,checked)
            checked = True
            html.append(u'<li><input %s checked/> ' % widgets.html_params(**options))
            html.append(u'<label for="%s">%s</label></li>' % (field_id, label))
        else:
            html.append(u'<tr><li><input %s /></tr> ' % widgets.html_params(**options))
            html.append(u'<tr><label for="%s">%s</label></li></tr>' % (field_id, label))
    html.append(u'</ul>')
    return u''.join(html)

def fetch_samples(id):
    cur = mysql.connection.cursor()
    cur.execute("SELECT * FROM samples WHERE id = %s", [id])
    return cur.fetchall()

class choose_samples(Form):
    samples_list = SelectMultipleField('Samples', choices=[], widget=select_multi_checkbox2)


#RUN PIPEPLINE
@app.route('/run_mosca/<string:id>', methods = ['GET', 'POST'])
@is_logged_in
def run_mosca(id):
    global state
    state = 'True'
    form = choose_samples(request.form)
    form.samples_list.choices = []
    samples = []
    cur = mysql.connection.cursor()
    project = cur.execute("SELECT * FROM projects WHERE id=%s",[id])
    project = cur.fetchone()
    n = project['author'].split(' ')
    tag_n= ''
    for i in n:
        tag_n += i[0]
    name = project['project_name']+'_'+tag_n+'_'+str(project['date']).split(' ')[0]
    #print(name)



    fetch = fetch_samples(id)
    if len(fetch) == 0:
        flash('Data needs to be added 1st' , 'danger')
        return redirect(url_for('dashboard'))
    samples_id = ''
    for row in fetch_samples(id):
        #print(row)
        samples.append(row)
        samples_id = '_'.join(map(str,str(row['s_id'])))

        #print(row['samples_name'])

        opt = 'Samples name: {} | Samples Condition: {} | Metagenomics Samples: {}; {} | Metatranscriptomics Samples {} {} '.format(str(row['samples_name']),str(row['samples_condition']),str(row['mg_sample1']),str(row['mg_sample2']),str(row['mt_sample1']),str(row['mt_sample2']))
        form.samples_list.choices.append((str(row),row['samples_name']))



        if request.method == 'POST' and form.validate():
            if not form.samples_list.data:
                flash('You need to select a sample' , 'danger')
                return redirect(request.url)
            samples2 = form.samples_list.data
            #print(samples2)
            #print(len(samples2))

            #execute

            for i in range(len(samples2)):

                cur.execute('INSERT INTO exe_projects(samples,id,s_id) VALUES(%s,%s,%s)',[samples2[i],int(samples[i]['id']),int(samples[i]['s_id'])])

            #commit to db
            mysql.connection.commit()
            #close connection
            cur.close()

            flash('Mosca is running' , 'success')
            name = name + '_' + samples_id
            exe_mosca=start_run(id,name,samples_id)

            #exe_mosca = start_run(id,name,samples_id)
            ### definir função para dar trigger no inicio da mosca
            return redirect(url_for('exe_mosca_pipe', id = id, name = name, samples_id=samples_id,exe_mosca=exe_mosca))

    cur.execute('DELETE FROM exe_projects')
    cur.execute('ALTER TABLE exe_projects AUTO_INCREMENT = 1')
        #pipe_samples = form.samples_list.choices
    #commit to db
    mysql.connection.commit()
    #close connection
    cur.close()
    return render_template('run_mosca.html', form=form, samples = samples)

#####################################################################################################################################################
#@app.route('/start',methods = ['GET', 'POST'])
def get_shell_script_output_using_communicate():
    #subprocess.Popen(['chmod','-x','execute_mosca2'])
    session = subprocess.Popen(['./execute_mosca2.sh'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = session.communicate()
    #print(stdout.decode('utf-8'))
    #if stderr:
        #print(str(Exception("Error "+str(stderr))))
        #raise Exception("Error "+str(stderr))

    return stdout.decode('utf-8')

def get_shell_script_output_using_check_output():
    stdout = check_output(['./execute_mosca2.sh']).decode('utf-8')
    return stdout

def create_project_directory(dir):
    subprocess.Popen(['mkdir',dir])

def subprocess_with_exp(expression):
    print(expression)
    print(expression.split('\t'))
    session = subprocess.Popen(expression.split('\t'), stdout=PIPE, stderr=PIPE)
    stdout, stderr = session.communicate()
    #print(stdout.decode('utf-8'))
    #if stderr:
        #return str(Exception("Error "+str(stderr)))
        #raise Exception("Error "+str(stderr))
    return stdout.decode('utf-8')



#PIPELINE EXECUTION + MONITORING
@app.route('/exe_mosca_pipe/<path:id>/<path:name>/<path:samples_id>/<path:exe_mosca>', methods = ['GET', 'POST'])
@is_logged_in
def exe_mosca_pipe(id, name, samples_id,exe_mosca):
    global state
    output = open('file.txt','r')
    steps = []
    report_out = []
    bin_ab = [[],[]]
    bin_sum = [[],[]]
    f_files = []
    if state == 'True':
        state = 'False'
        print('\n'+state+'\n')

        return render_template('exe_mosca_pipe.html', id=id, name = name, samples_id=samples_id, exe_mosca=exe_mosca, steps=steps, report_out=report_out, f_files=f_files, bin_ab=bin_ab,bin_sum=bin_sum,state=state)

        print('\n'+state+'\n')

    elif state == 'False':
        state = 'Pass'
        print('\n'+state+'\n')

        #subprocess.run(exe_mosca.split('\t'), stdout=PIPE, check = True)
        #get_shell_script_output_using_communicate()
        subprocess.run(exe_mosca.split('\t'),stdout=PIPE,stderr=subprocess.STDOUT,capture_output=True)

    else:
        print('\n'+state+'\n')

        for line in output:
            print(line)
            steps.append(line.rstrip('\n'))
            if 'assembly' in line :
                report = open('static/{}/Assembly/quality_control/report.tsv'.format(name), 'r')
                for l in report:
                    report_out.append(l.rstrip('\n'))




            if 'binning' in line:
                #bin_files = get_filelist(os.path.join(app.instance_path)[0:-8]+'static/Binning')

                bin_files = glob.glob(os.path.join(app.instance_path)[0:-8]+'static/{}/Binning/markerset*'.format(name))
                #static/Binning/markerset40.summary'abundance
                print(bin_files)
                for i in range(len(bin_files)):
                    if bin_files[i] != 'None' :
                        if 'abundance' in bin_files[i] :
                            f = open('static'+ bin_files[i].split('static')[1],'r')
                            #f = pd.read_csv(bin_files[i][0][34])
                            print(f)

                            bin_ab[0].append(bin_files[i][-36:].split('/')[-1])
                            for l in f:

                                bin_ab[1].append(l.rstrip('\n'))
                        else:
                            f = open('static'+ bin_files[i].split('static')[1],'r')
                            bin_sum[0].append(bin_files[i][-34:].split('/')[-1])
                            for l in f:

                                bin_sum[1].append(l.rstrip('\n'))
                        #f = open()
                print(bin_ab)
                print(bin_sum)






            if 'preprocessing' in line:
                #pre_files = get_filelist(os.path.join(app.instance_path)[0:-8]+'static/Preprocess/FastQC')
                pre_files = glob.glob(os.path.join(app.instance_path)[0:-8]+'static/{}/Preprocess/FastQC/quality_trimmed_*_paired_fastqc.html'.format(name))

                for file in pre_files:
                    #print(file.split('static')[1])

                    f_files.append('/static'+file.split('static')[1])


                print('hi',pre_files)
                print(f_files)
                state = 'True'

        print('HI\n',report_out)
        print(os.path.join(app.instance_path)[0:-8])

        return render_template('exe_mosca_pipe.html', id=id, name = name, samples_id=samples_id, exe_mosca=exe_mosca, steps=steps, report_out=report_out, f_files=f_files, bin_ab=bin_ab,bin_sum=bin_sum,state=state)

###############################################################
@app.route('/annotation/taxonomy/<path:name>', methods = ['GET', 'POST'])
def krona1(name):
    return send_file('static/{}/Annotation/mg_taxonomy.html'.format(name))

@app.route('/annotation/cogs/<path:name>', methods = ['GET', 'POST'])
def krona2(name):
    return send_file('static/{}/Annotation/mg_cogs.html'.format(name))


#EXPRESSION CONSTRUCTION
def start_run(id, name, samples_id):
    global state
    state = 'True'


    create_project_directory('static/{}'.format(name))

    cur = mysql.connection.cursor()
    project = cur.execute("SELECT * FROM projects WHERE id=%s",[id])
    project = cur.fetchone()
    samples = cur.execute("SELECT samples FROM exe_projects WHERE id=%s",[id])
    samples = cur.fetchall()
    file_exp= ''


    for i in range(len(samples)):
        sample = ast.literal_eval(samples[i]['samples'])
        #print(sample)
        if i >=1:
            file_exp += ' '
        if sample['mg_sample1'] and sample['mg_sample2'] and sample['mt_sample1'] and sample['mt_sample2']:
            file_exp += sample['mg_sample1'] + ',' + sample['mg_sample2'] + ':' + sample['mt_sample1'] + ',' + sample['mt_sample2']
        elif sample['mg_sample1'] and sample['mg_sample2'] and sample['mt_sample1'] and not sample['mt_sample2']:
            file_exp += sample['mg_sample1'] + ',' + sample['mg_sample2'] + ':' + sample['mt_sample1']
        elif sample['mg_sample1'] and sample['mg_sample2'] and not sample['mt_sample1'] and not sample['mt_sample2']:
            file_exp += sample['mg_sample1'] + ',' + sample['mg_sample2']
        elif not sample['mg_sample1'] and not sample['mg_sample2'] and sample['mt_sample1'] and sample['mt_sample2']:
            file_exp += sample['mt_sample1'] + ',' + sample['mt_sample2']
        elif sample['mg_sample1'] and not sample['mg_sample2'] and sample['mt_sample1'] and not sample['mt_sample2']:
            file_exp += sample['mg_sample1'] + ':' + sample['mt_sample1']
        elif sample['mg_sample1'] and not sample['mg_sample2'] and not sample['mt_sample1'] and not sample['mt_sample2']:
            file_exp += sample['mg_sample1']
        elif not sample['mg_sample1'] and not sample['mg_sample2'] and  sample['mt_sample1'] and not sample['mt_sample2']:
            file_exp += sample['mt_sample1']

    exe_try = 'python MOSCA/scripts/mosca.py --files\t{}\t'.format(file_exp)
    exe_mosca = 'python MOSCA/scripts/mosca.py --files\t{}\t'.format(file_exp)

    print(project)
    if project['sequencing']=='PE':
        st_exp = '--sequencing-technology\tpaired'
    if project['sequencing']=='SE':
        st_exp = '--sequencing-technology\tsingle'

    ass_exp = '--assembler\t{}'.format(project['assembler'].lower())

    db_exp = '--annotation-database\t{}'.format(project['database_dir'])

    out_exp = '--output\t{}static/{}'.format(os.path.join(app.instance_path)[0:-8],name)

    no_exp = ''
    if not project['preprocessing']:
        no_exp += '\t--no-preprocessing'
    if not project['assembly']:
        no_exp += '\t--no-assembly'
    if not project['binning']:
        no_exp += '\t--no-binning'

    outlvl_exp = '--output-level\t{}'.format(project['output_lvl'])

    tod_exp = '--type-of-data\t{}'.format(project['data_type'])
    #if sample['mg_sample1'] and not sample['mt_sample1']:
    #    tod_exp = '--type-of-data\tmetagenomics'
    #if sample['mt_sample1'] and not sample['mg_sample1']:
    #    tod_exp = '--type-of-data\tmetatranscriptomics'
    #if sample['mt_sample1'] and sample['mg_sample1']:
    #    tod_exp = '--type-of-data\tmetagenomics,metatranscriptomics'

    dict_s = []
    for i in samples:
        dict_s.append(literal_eval(i['samples']))
    #print(dict_s)

    res = []
    for s in dict_s:
        res.append(s['samples_condition'])
    scond_exp=','.join(res)
    #print(scond_exp)
    scond_exp='--conditions\t'+scond_exp

    thr_exp = '--threads\t{}'.format(project['threads'])

    memory_exp = '--memory\t{}'.format(project['memory'])

    gene_set_exp = '--marker-gene-set\t{}'.format(project['marker'].rstrip(' '))

    #exe_mosca = 'python MOSCA/scripts/mosca.py\t--files\t{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(file_exp,st_exp,ass_exp,db_exp,out_exp,no_exp,outlvl_exp,tod_exp,scond_exp,thr_exp,memory_exp,gene_set_exp)
    mosca_exe = 'MOSCA/scripts/mosca.py'
    #mosca_exe = '/mnt/HDDStorage/jsequeira/MOSCA/scripts/mosca.py'

    exe_mosca = 'python\t{}\t--files\t{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(mosca_exe,file_exp,st_exp,ass_exp,db_exp,out_exp,no_exp,outlvl_exp,tod_exp,scond_exp,thr_exp,memory_exp,gene_set_exp)
    #(exe_mosca.split(subprocess.run'\t'), stdout=PIPE, check = True)

    #exe_mosca = 'python MOSCA/scripts/mosca.py --files {} --output-dir output_directory'.format(file_exp)
    file = open('execute_mosca2.sh','w')
    file.write('#!/bin/sh'+'\n')
    #file.write('chmod -x execute_mosca2.sh'+'\n')
    #ile.write('echo ' + "'{}'".format(exe_mosca))
    file.write("{}".format(exe_mosca))
    #file.write("{}".format(exe_try))
    file.close()


    #print(exe_mosca)

    #exe_mosca = get_shell_script_output_using_check_output()
    #exe_mosca = get_shell_script_output_using_communicate()
    #print(exe_mosca)
    #print(background_process())
    #get_shell_script_output_using_communicate()
    #print(os.path.join(app.instance_path))
    #get_shell_script_output_using_communicate()
    #subprocess_with_exp(exe_mosca)

    return exe_mosca



@app.route('/background_process')
def background_process():
    while True:
        try:
            with open('file.txt', 'r') as file:
                step = 0
                steps = {}
                msgs = []
                if os.stat("file.txt").st_size == 0:
                    steps[0]='Mosca is starting'
                else:
                    steps[0]='Mosca is starting'
                    for line in file:
                        step += 1
                        steps[step]=line
                if steps == 5:
                    msgs.append(steps)
                    break

                else:
                    msgs.append(steps)
                    #return msgs
            return msgs
        except Exception as e:
            raise



if __name__ == '__main__':
    app.run(host='0.0.0.0')

#def main():
    #app.run(host='0.0.0.0')

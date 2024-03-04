scATAC Walkthrough
===================

About
-------------------------
ABC can predict cell-type-specific enhancer-gene interactions using only the chromatin accessibility profile obtained through scATAC-seq. Given the sparsity of scATAC data, we recommend pseudo-bulking the scATAC fragments by cell cluster. Each pseudo-bulked scATAC profile should contain at least 3 million (preferably 6 million) unique fragments to ensure optimal performance.

We'll walk through an example of running a scATAC dataset through ABC. 

1. Converting to TagAlign
-------------------------
.. admonition:: **Prerequisite**

    Make sure you've followed the environment setup steps in :ref:`GettingStarted` and activate your abc conda environment.

When working with scATAC files, you may start off with a `fragment file <https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/fragments?src=social&lss=facebook&cnm=soc-fb-ra_g-program-fb-ra_g-program&cid=7011P000000y072>`_. Since ABC works with `tagAlign file format <https://genome.ucsc.edu/FAQ/FAQformat.html#format15>`_,  we'll want to convert our fragment file to tagAlign. If you have a tagAlign file already, great! You can skip to the stage where we create the biosamples table.

Let's download an example fragment file from the internet

.. code-block:: bash

    (abc-env) [atan5133@sh03-04n23 /oak/stanford/groups/engreitz/Users/atan5133/data] (job 37050919) $ wget https://www.encodeproject.org/files/ENCFF794UXO/@@download/ENCFF794UXO.tar.gz
    
    (abc-env) [atan5133@sh03-04n23 /oak/stanford/groups/engreitz/Users/atan5133/data] (job 37050919) $ tar -xf ENCFF794UXO.tar.gz

    (abc-env) [atan5133@sh03-04n23 /oak/stanford/groups/engreitz/Users/atan5133/data] (job 37050919) $ ls encode_scatac_dcc_2/results/ENCSR308ZGJ-1/fragments/
        fragments.tsv.gz  fragments.tsv.gz.tbi


Let's convert the fragment file to .tagAlign and index it 

.. code-block:: bash

    (abc-env) [atan5133@sh03-04n23 /oak/stanford/groups/engreitz/Users/atan5133/data] (job 37050919) $ cd encode_scatac_dcc_2/results/ENCSR308ZGJ-1/fragments
    
    (abc-env) [atan5133@sh03-04n23 /oak/stanford/groups/engreitz/Users/atan5133/data/encode_scatac_dcc_2/results/ENCSR308ZGJ-1/fragments] (job 37050919) $ LC_ALL=C zcat fragments.tsv.gz | sed '/^#/d' | awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"}' | sort -k 1,1V -k 2,2n -k3,3n --parallel 5 | bgzip -c > tagAlign.gz  # Adjust --parallel 5 based on number of cpus you have. The more cpus, the faster

    (abc-env) [atan5133@sh03-04n24 /oak/stanford/groups/engreitz/Users/atan5133/data/encode_scatac_dcc_2/results/ENCSR308ZGJ-1/fragments] (job 37151429) $ tabix -p bed tagAlign.gz

Sorting the tagAlign typically takes a while and varies based on the size of the file. Setting LC_ALL=C and adding parallelization makes the sorting much faster. 
    

2. Creating Your Config
-------------------------

Once we have the tagAlign file, we'll create our biosample_table config. Let's assume we don't have any HiC data, so we'll opt to use powerlaw as our estimate for contact.
Our biosamples config would look like the following:

.. list-table::
   :header-rows: 1
   :widths: auto

   * - biosample
     - DHS
     - ATAC
     - H3K27ac
     - default_accessibility_feature
     - HiC_file
     - HiC_type
     - HiC_resolution
     - alt_TSS
     - alt_genes
   * - K562
     - 
     - /oak/stanford/groups/engreitz/Users/atan5133/data/encode_scatac_dcc_2/results/ENCSR308ZGJ-1/fragments/tagAlign.gz
     - 
     - ATAC
     - 
     -
     -
     - 
     - 

Update the `ABC Config <https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/dev/config/config.yaml#L5>`_ to point the biosample tsv file you just created


3. Running ABC
---------------

.. code-block:: bash

    (abc-env) [atan5133@sh02-ln01 login /oak/stanford/groups/engreitz/Users/atan5133/
    ABC-Enhancer-Gene-Prediction]$ snakemake -j1


All your results will go to the ``results`` directory of your ABC directory! 
The actual predictions will be stored under the ``results/{biosample}/Predictions`` folder.

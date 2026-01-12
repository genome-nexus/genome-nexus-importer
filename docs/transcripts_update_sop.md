Check list for transcript updating:
-	Work with OncoKB and clinical team to get an approval on the change
    -	The new transcript needs to be on Ensembl111 and has sub version
    -	Genome Nexus and OncoKB should make the updates at the same time.
-	Update Genome Nexus Importer mskcc isoform list, document changes in the changelog
-	Making a new Genome Nexus mongodb release, update all Genome Nexus instances to use the latest image
    -	annotation.genomenexus.org
    -	www.genomenexus.org
-	Find related samples, get samples ids and work with pipeline team to requeue the samples.
-	Go to cbioportal to check if the default transcript is switched to new transcript, and if requeued samples are on new annotations.

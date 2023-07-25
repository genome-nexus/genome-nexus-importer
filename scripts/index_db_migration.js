// Command to run script: "mongo < index_db_migration.js"
db.adminCommand('listDatabases');
db = db.getSiblingDB('annotator');
db.getCollectionNames();
searchCollection = db.getCollection('index');
const effectPriority = {};
effectPriority["transcript_ablation"] = 1; // A feature ablation whereby the deleted region includes a transcript feature
effectPriority["exon_loss_variant"] = 1; // A sequence variant whereby an exon is lost from the transcript
effectPriority["splice_donor_variant"] = 2; // A splice variant that changes the 2 base region at the 5" end of an intron
effectPriority["splice_acceptor_variant"] = 2; // A splice variant that changes the 2 base region at the 3" end of an intron
effectPriority["stop_gained"] = 3; // A sequence variant whereby at least one base of a codon is changed] = resulting in a premature stop codon, leading to a shortened transcript
effectPriority["frameshift_variant"] = 3; // A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
effectPriority["stop_lost"] = 3; // A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
effectPriority["start_lost"] = 4; // A codon variant that changes at least one base of the canonical start codon
effectPriority["initiator_codon_variant"] = 4; // A codon variant that changes at least one base of the first codon of a transcript
effectPriority["disruptive_inframe_insertion"] = 5; // An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
effectPriority["disruptive_inframe_deletion"] = 5; // An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
effectPriority["inframe_insertion"] = 5; // An inframe non synonymous variant that inserts bases into the coding sequence
effectPriority["inframe_deletion"] = 5; // An inframe non synonymous variant that deletes bases from the coding sequence
effectPriority["protein_altering_variant"] = 5; // A sequence variant which is predicted to change the protein encoded in the coding sequence
effectPriority["missense_variant"] = 6; // A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
effectPriority["conservative_missense_variant"] = 6; // A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
effectPriority["rare_amino_acid_variant"] = 6; // A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
effectPriority["transcript_amplification"] = 7; // A feature amplification of a region containing a transcript
effectPriority["splice_region_variant"] = 8; // A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
effectPriority["start_retained_variant"] = 9; // A sequence variant where at least one base in the start codon is changed, but the start remains
effectPriority["stop_retained_variant"] = 9; // A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
effectPriority["synonymous_variant"] = 9; // A sequence variant where there is no resulting change to the encoded amino acid
effectPriority["incomplete_terminal_codon_variant"] = 10; // A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
effectPriority["coding_sequence_variant"] = 11; // A sequence variant that changes the coding sequence
effectPriority["mature_mirna_variant"] = 11; // A transcript variant located with the sequence of the mature miRNA
effectPriority["exon_variant"] = 11; // A sequence variant that changes exon sequence
effectPriority["5_prime_utr_variant"] = 12; // A UTR variant of the 5" UTR
effectPriority["5_prime_utr_premature_start_codon_gain_variant"] = 12; // snpEff-specific effect, creating a start codon in 5" UTR
effectPriority["3_prime_utr_variant"] = 12; // A UTR variant of the 3" UTR
effectPriority["non_coding_exon_variant"] = 13; // A sequence variant that changes non-coding exon sequence
effectPriority["non_coding_transcript_exon_variant"] = 13; // snpEff-specific synonym for non_coding_exon_variant
effectPriority["non_coding_transcript_variant"] = 14; // A transcript variant of a non coding RNA gene
effectPriority["nc_transcript_variant"] = 14; // A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
effectPriority["intron_variant"] = 14; // A transcript variant occurring within an intron
effectPriority["intragenic_variant"] = 14; // A variant that occurs within a gene but falls outside of all transcript features. This occurs when alternate transcripts of a gene do not share overlapping sequence
effectPriority["intragenic"] = 14; // snpEff-specific synonym of intragenic_variant
effectPriority["nmd_transcript_variant"] = 15; // A variant in a transcript that is the target of NMD
effectPriority["upstream_gene_variant"] = 16; // A sequence variant located 5" of a gene
effectPriority["downstream_gene_variant"] = 16; // A sequence variant located 3" of a gene
effectPriority["tfbs_ablation"] = 17; // A feature ablation whereby the deleted region includes a transcription factor binding site
effectPriority["tfbs_amplification"] = 17; // A feature amplification of a region containing a transcription factor binding site
effectPriority["tf_binding_site_variant"] = 17; // A sequence variant located within a transcription factor binding site
effectPriority["regulatory_region_ablation"] = 17; // A feature ablation whereby the deleted region includes a regulatory region
effectPriority["regulatory_region_amplification"] = 17; // A feature amplification of a region containing a regulatory region
effectPriority["regulatory_region_variant"] = 17; // A sequence variant located within a regulatory region
effectPriority["regulatory_region"] = 17; // snpEff-specific effect that should really be regulatory_region_variant
effectPriority["feature_elongation"] = 18; // A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
effectPriority["feature_truncation"] = 18; // A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
effectPriority["intergenic_variant"] = 19; // A sequence variant located in the intergenic region, between genes
effectPriority["intergenic_region"] = 19; // snpEff-specific effect that should really be intergenic_variant
effectPriority[""] = 20;
const AA3TO1 = [
    ['Ala', 'A'], ['Arg', 'R'], ['Asn', 'N'], ['Asp', 'D'], ['Asx', 'B'], ['Cys', 'C'],
    ['Glu', 'E'], ['Gln', 'Q'], ['Glx', 'Z'], ['Gly', 'G'], ['His', 'H'], ['Ile', 'I'],
    ['Leu', 'L'], ['Lys', 'K'], ['Met', 'M'], ['Phe', 'F'], ['Pro', 'P'], ['Ser', 'S'],
    ['Thr', 'T'], ['Trp', 'W'], ['Tyr', 'Y'], ['Val', 'V'], ['Xxx', 'X'], ['Ter', '*']
];
db.getCollection('vep.annotation').find().forEach( function(annotation) {
    let record = {
        _id: annotation._id,
        variant: annotation._id,
        hugoSymbol: [],
        hgvsp: [],
        hgvsc: [],
        cdna: [],
        hgvspShort: [],
        rsid: [],
    };
    if (annotation.transcript_consequences && annotation.transcript_consequences.length > 0) {
        for (let transcriptConsequence of annotation.transcript_consequences) {
            // Hugo symbol
            if (transcriptConsequence.gene_symbol && !record.hugoSymbol.includes(transcriptConsequence.gene_symbol)) {
                record.hugoSymbol.push(transcriptConsequence.gene_symbol);
            }
            // hgvsp
            if (transcriptConsequence.hgvsp && !record.hgvsp.includes(transcriptConsequence.hgvsp)) {
                record.hgvsp.push(transcriptConsequence.hgvsp);
            }
            // hgvsc
            if (transcriptConsequence.hgvsc && !record.hgvsc.includes(transcriptConsequence.hgvsc)) {
                record.hgvsc.push(transcriptConsequence.hgvsc);
            }
            // cdna
            if (transcriptConsequence.hgvsc && transcriptConsequence.hgvsc.includes(':c.')) {
                const cdna = transcriptConsequence.hgvsc.split(':')[1];
                if (!record.cdna.includes(cdna)) {
                    record.cdna.push(cdna);
                }
            }
            // hgvspShort
            if (transcriptConsequence.hgvsp) {
                let hgvspShort = null;
                if (transcriptConsequence.hgvsp.includes(':p.')) {
                    let highestPriorityConsequence = null;
                    let highestPriority = 21; // highest number in effectPriority map is 20
                    // check if this variant is splice
                    if (transcriptConsequence.consequenceTerms) {
                        transcriptConsequence.consequenceTerms.forEach((consequenceTerm) => {
                                if ((effectPriority[consequenceTerm] || highestPriority) < highestPriority) {
                                    highestPriorityConsequence = consequenceTerm;
                                    highestPriority = effectPriority[consequenceTerm] || highestPriority;
                                }
                            }
                        );
                    }
                    let variantClassification = highestPriorityConsequence || annotation.mostSevereConsequence;
                    // if variant is not splice, resolve from hgvsp
                    if (!(variantClassification && variantClassification.toLowerCase().contains("splice"))) {
                        hgvspShort = transcriptConsequence.hgvsp.split(':')[1];
                        for (let i = 0; i < 24; i++) {
                            if (hgvspShort.includes(AA3TO1[i][0])) {
                                let replace = AA3TO1[i][0];
                                let re = new RegExp(replace, "g");
                                hgvspShort = hgvspShort.replace(re, AA3TO1[i][1]);
                            }
                        }
                    }             
                }
                // skip other cases
                if (!record.hgvspShort.includes(hgvspShort) && hgvspShort!= null) {
                    record.hgvspShort.push(hgvspShort);
                }
            }
        }
    }
    searchCollection.insert(record);
})
printjson("Migration done!");
db.getCollection('index').createIndex({'variant':1});
db.getCollection('index').createIndex({'hugoSymbol':1});
db.getCollection('index').createIndex({'hgvspShort':1});
db.getCollection('index').createIndex({'hgvsp':1});
db.getCollection('index').createIndex({'cdna':1});
db.getCollection('index').createIndex({'hgvsc':1});
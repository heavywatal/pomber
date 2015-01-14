library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pipeR)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))

library(Biostrings)
library(GenomicRanges)

.cache_dir = '~/Dropbox/cache'
setwd('~/git/pomber')

make_windows = function(len, width, step)
    successiveIRanges(rep(width, (len - width) / step + 1), step - width)

#########1#########2#########3#########4#########5#########6#########7#########
## polymorphism oligo_stats

.cache_file = file.path(.cache_dir, 'origo_stats.rda')
.cache_file = 'origo_stats.rda'
if (file.exists(.cache_file)) {
    load(.cache_file)
} else {
    .data_dir = '~/db/ilabo/pombe'
    .patt = '.*([^/]+)\\.bin\\.fa\\.gz$'
    .seq_files = list.files(.data_dir, .patt, full.names=TRUE) %>>% (? head(.))
    names(.seq_files) = str_match(.seq_files, .patt)[,2]
    .bssl = llply(.seq_files, readDNAStringSet) %>>% DNAStringSetList()

    ldply(.bssl %>>% head(), alphabetFrequency) %>>% (? head(.)) %>>%
        filter(rep(c(TRUE, FALSE), length(.) / 2)) %>>%
        summarise(AT=sum(A+T), GC=sum(G+C)) %>>%
        mutate(ratio=AT / (AT+GC))
    # 0.6284053

    .patt = '.*?([^/]+)\\.pi-(\\d+)-1\\.txt\\.gz$'
    .stat_files = list.files(.data_dir, .patt, full.names=TRUE) %>>% (? head(.))
    read_stats = function(.infile) {
        message(.infile)
        .match = str_match(.infile, .patt)
        chromosome = .match[2]
        .width = as.integer(.match[3])
        .swstats = read.delim(.infile)
        if (.width > 1) {
            .swstats = .swstats %>>% head(1 - .width)
        }
        .bss = .bssl[[chromosome]]
        .windows = make_windows(length(.bss[[1]]), .width, 1)
        data.frame(
            major=Views(.bss[[1]], .windows) %>>% as.character(),
            minor=Views(.bss[[2]], .windows) %>>% as.character(),
            stringsAsFactors=FALSE) %>>%
            cbind(.swstats) %>>%
            filter(!str_detect(major, 'N') & !str_detect(minor, 'N'))
    }
    .raw = ldply(.stat_files, read_stats, .parallel=TRUE)
    .raw %>>% sample_n(30)

    oligo_stats = .raw
    if (FALSE) {  # divergence, not polymorphism
        oligo_stats = oligo_stats %>>% dplyr::select(-theta_pi, -theta_w, -TajimasD)
    }
    save(oligo_stats, file=.cache_file)
}

#########1#########2#########3#########4#########5#########6#########7#########
## polymorphism table

.outfile = 'polymorphism/pombe-mean.tsv.gz'

.major = oligo_stats %>>% tbl_df() %>>%
    dplyr::select(-minor, -theta_w) %>>%
    plyr::rename(c(major='oligo')) %>>%
    mutate(oligo=oligo %>>% DNAStringSet() %>>% reverseComplement() %>>% as.character() %>>% pmin(oligo))

.freq = .major %>>% group_by(oligo) %>>% tally()

.mean = .major %>>%
    group_by(oligo) %>>%
    summarise_each(funs(mean(., na.rm=TRUE))) %>>%
    left_join(x=.freq, by='oligo') %>>%
    mutate(oligo=paste(oligo,
        oligo %>>% DNAStringSet() %>>% reverseComplement() %>>% as.character(),
        sep='/'))

write.table(.mean, gzfile(.outfile), sep='\t', row.names=FALSE)

.mean %>>%
    group_by(width=nchar(oligo)) %>>%
    arrange(width, theta_pi)

.mean %>>%
    group_by(width=nchar(oligo)) %>>%
    summarise_each(funs(quantile(., probs=0.10)), -oligo)

#########1#########2#########3#########4#########5#########6#########7#########
## evenness

.outfile = 'evenness/coefvar-ref.csv.gz'

chromosomes = list.files('distribution',
                         pattern='.*\\.chromosome\\.I+\\.fa\\.gz',
                         full.names=TRUE) %>>%
              readDNAStringSet()
names(chromosomes) = names(chromosomes) %>>% str_extract('^I+')

revcomp = function(x)
    x %>>% DNAStringSet() %>>% reverseComplement() %>>% as.character()

complement = function(x)
    chartr('ATGC', 'TACG', x)

make_query = function(s) {
    rc = revcomp(s)
    if (s == rc) {
        PDict(s)
    } else {
        PDict(c(s, rc))
    }
}

.weight = .weight %>>% select(matches('[ATGC]'))
.null = paste(sample(names(.weight), 2000000, replace=TRUE, .weight), collapse='')
chromosomes = c(chromosomes, DNAStringSet(c(randomized=.null)))

sliding_window = function(.query, .width, .step) {
    .data = ldply(chromosomes, function(.chr) {
        .windows = make_windows(length(.chr), .width, .step)
        .hit = matchPDict(.query, .chr) %>>% Biostrings::unlist()
        data.frame(
            observed=countOverlaps(.windows, .hit),
            stringsAsFactors=FALSE
        )
    }, .id='chr')
}

if (FALSE) {
    .query = 'ATGCGC'
    .data = sliding_window(make_query(.query), .width, .step)
    .data %>>%
        group_by(chr) %>>%
        summarise(mean=mean(observed), var=var(ovserved), sd=sd(observed), vm_ratio=var/mean, cv=sd/mean)
}

.oligos = strsplit('ATGC', '') %>>%
    rep(6) %>>%
    expand.grid() %>>%
    mlply(function(...) paste0(...)) %>>%
    unlist(use.names=FALSE)
names(.oligos) = .oligos

.queries = data_frame(motif=.oligos, revcomp=revcomp(.oligos)) %>>%
    mutate(motif=pmin(motif, revcomp)) %>>%
    select(-revcomp) %>>%
    filter(!duplicated(motif)) %>>%
    mutate(revcomp=revcomp(motif), motif=ifelse(motif==revcomp, motif, paste(motif, revcomp, sep='|'))) %>>%
    select(-revcomp)

.width = 20000
.step = 20000

.raw = mdply(.queries, function(motif) {
    .query = strsplit(motif, '\\|') %>>% unlist() %>>% PDict()
    sliding_window(.query, .width, .step)
}, .id='motif', .parallel=TRUE)

.each_chr = .raw %>>%
    group_by(motif, chr) %>>%
    summarise(mean=mean(observed), var=var(observed), sd=sd(observed), vm_ratio=var/mean, cv=sd/mean)

.total_genome = .raw %>>%
    filter(grepl('I', chr)) %>>%
    group_by(motif) %>>%
    summarise(mean=mean(observed), var=var(observed), sd=sd(observed), vm_ratio=var/mean, cv=sd/mean) %>>%
    mutate(chr='total')

.data = bind_rows(.each_chr, .total_genome) %>>%
    arrange(motif)

if (FALSE) {
    .data %>>%
        filter(grepl('I', chr)) %>>%
        arrange(vm_ratio) %>%
        each(head, tail)(.)
}

write.table(.data, gzfile(.outfile), sep='\t', row.names=FALSE, quote=FALSE)

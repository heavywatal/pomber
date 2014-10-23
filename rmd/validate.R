library(Biostrings)
library(GenomicRanges)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))

revcomp = function(x)
    x %>>% DNAStringSet() %>>% reverseComplement() %>>% as.character()

load('data/chromosomes.rda')
#########1#########2#########3#########4#########5#########6#########7#########

.chr = 1
.start = 200000
.width = 100000
.remove_N = FALSE
.shuffle = FALSE

.seq = chromosomes[[.chr]][.start:(.start + .width - 1)]

if (.remove_N) {
    .chars = .seq %>>% as.character() %>>% strsplit('') %>>% unlist()
    .seq = .chars[.chars != 'N'] %>>% paste(collapse='') %>>% DNAString() %>>% (? .)
}
if (.shuffle) {
    .seq = .seq[sample(length(.seq))] %>>% (? .)
}

.nucs = c('A', 'C', 'G', 'T', 'N')
.nucl_freq = .seq %>>%
    alphabetFrequency(as.prob=TRUE) %>>%
    `[`(.nucs) %>>% (? .)
.motifs = expand.grid(.nucs %>>% list() %>>% rep(4)) %>>%
    mlply(paste0) %>>% unlist() %>>% (? head(.))
names(.motifs) = .motifs

.data = ldply(.motifs, function(.query) {
    rc = revcomp(.query)
    .obs = if (.query == rc) {
        .q = .query
        countPattern(.query, .seq)
    } else {
        .q = c(.query, rc)
        countPattern(.query, .seq) +
        countPattern(rc, .seq)
    }
    .prob = .q %>>%
        strsplit('') %>>%
        llply(paste, collapse='*') %>>%
        unlist() %>>%
        paste(collapse='+')# %>>% (? .)
    .gc = sum(alphabetFrequency(DNAString(.query))[c('G', 'C')])
    .nucl_freq %>>% t() %>>% data.frame() %>>%
        transmute_(p=.prob) %>>%
        mutate(expected = p * (.width - nchar(.query) + 1)) %>>%
        mutate(observed = .obs, GC=.gc)
}, .id='query', .parallel=TRUE) %>>% tbl_df() %>>% (? .)

.data %>>% summarise_each(funs(sum))
.data %>>% filter(observed < expected)

.range = range(.data$expected) %>>% (? .)
.data %>>%
    ggplot(aes(expected, observed))+
    geom_point(aes(colour=GC), alpha=0.4)+
    geom_line(data=data_frame(expected=.range, observed=.range), colour='#FF3300')+
    theme_bw()


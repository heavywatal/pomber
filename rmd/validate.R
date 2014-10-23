library(Biostrings)
library(GenomicRanges)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))

revcomp = function(x)
    x %>>% DNAStringSet() %>>% reverseComplement() %>>% as.character()

make_query = function(s) {
    rc = revcomp(s)
    if (s == rc) {
        PDict(s)
    } else {
        PDict(c(s, rc))
    }
}

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

.nucs = c('A', 'C', 'G', 'T')
.nucl_freq = .seq %>>%
    alphabetFrequency(as.prob=TRUE) %>>%
    `[`(.nucs) %>>% (? .)
.motifs = expand.grid(.nucs %>>% list() %>>% rep(4)) %>>%
    mlply(paste0) %>>% unlist() %>>% (? head(.))
names(.motifs) = .motifs

.data = ldply(.motifs, function(s) {
    .query = make_query(s)
    .prob = .query@dict0 %>>%
        as.character() %>>%
        strsplit('') %>>%
        llply(paste, collapse='*') %>>%
        unlist() %>>%
        paste(collapse='+')# %>>% (? .)
    .gc = sum(alphabetFrequency(DNAString(s))[c('G', 'C')])
    .nucl_freq %>>% t() %>>% data.frame() %>>%
        transmute_(p=.prob) %>>%
        mutate(expected = p * (.width - width(.query)[1] + 1)) %>>%
        mutate(observed = sum(countPDict(.query, .seq)), GC=.gc)
}, .id='query', .parallel=TRUE) %>>% tbl_df() %>>% (? .)

.data %>>% summarise_each(funs(sum))

.range = range(.data$expected) %>>% (? .)
.data %>>%
    ggplot(aes(expected, observed))+
    geom_point(aes(colour=GC), alpha=0.4)+
    geom_line(data=data_frame(expected=.range, observed=.range), colour='#FF3300')+
    theme_bw()


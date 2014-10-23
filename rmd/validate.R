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
.seq %>>% as.character()

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

.raw = ldply(.motifs, function(.query) {
    .obs = countPattern(.query, .seq)
    .prob = .query %>>%
        strsplit('') %>>%
        llply(paste, collapse='*') %>>%
        unlist() %>>%
        paste(collapse='+')# %>>% (? .)
    .qalpha = alphabetFrequency(DNAString(.query))
    .nucl_freq %>>% t() %>>% data.frame() %>>%
        transmute_(p=.prob) %>>%
        mutate(expected = p * (.width - nchar(.query) + 1)) %>>%
        mutate(observed = .obs, GC=sum(.qalpha[c('G', 'C')]), N=sum(.qalpha['N']))
}, .id='query', .parallel=TRUE) %>>% tbl_df() %>>% (? .)

.data = .raw %>>%
    mutate(query=as.character(query)) %>>%
    arrange(query) %>>%
    mutate(query=query %>>% revcomp() %>>% pmin(query)) %>>%
    group_by(query) %>>%
    summarise_each(funs(sum)) %>>%
    mutate(GC=GC/2, N=N/2) %>>%
    (? summarise_each(., funs(sum), -query))

.data %>>% arrange(desc(observed))
.data %>>% filter(observed < expected)
.data %>>% filter(observed > expected)

.range = range(.rc$expected) %>>% (? .)
.data %>>%
    ggplot(aes(expected, observed))+
    geom_point(aes(colour=N), alpha=0.4)+
    geom_line(data=data_frame(expected=.range, observed=.range), colour='#FF3300')+
    theme_bw()

library(Biostrings)
library(GenomicRanges)
library(doMC)
doMC::registerDoMC(min(parallel::detectCores(), 12))

revcomp = function(x)
    x %>>% DNAStringSet() %>>% reverseComplement() %>>% as.character()

load('data/ensembl.rda')
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
.motifs = expand.grid(.nucs %>>% list() %>>% rep(5)) %>>%
    mlply(paste0) %>>% unlist() %>>% (? head(.))
names(.motifs) = .motifs

.raw = ldply(.motifs, function(.query) {
    .obs = countPattern(.query, .seq)
    .formula = strsplit(.query, '') %>>%
        unlist() %>>%
        paste(collapse='*')
    .nucl_freq %>>% t() %>>% data.frame() %>>%
        transmute_(p=.formula) %>>%
        mutate(expected = p * (.width - nchar(.query) + 1),
               observed = .obs)
}, .id='query', .parallel=TRUE) %>>% tbl_df() %>>% (? .)

.data = .raw %>>%
    mutate(query=as.character(query)) %>>%
    arrange(query) %>>%
    mutate(query=pmin(query, query %>>% revcomp())) %>>%
    group_by(query) %>>%
    summarise_each(funs(sum)) %>>%
    cbind(alphabetFrequency(DNAStringSet(.$query)) %>>%
        data.frame() %>>%
        select(matches('[ACGTN]'))) %>>%
    tbl_df() %>>%
    (? summarise_each(., funs(sum), -query))

.data %>>% arrange(desc(observed))

.data %>>%
    group_by(XXXXN=N>0) %>>%
    tally(wt=p)
(1.0 - .nucl_freq['N']) ^ 5

.range = range(.data$expected) %>>% (? .)
.data %>>%
    filter(N<5) %>>%
    ggplot(aes(expected, observed))+
    geom_point(aes(colour=N), alpha=0.4)+
    geom_line(data=data_frame(expected=.range, observed=.range), colour='#FF3300')+
    theme_bw()+
    scale_colour_gradient(low='#00CCCC', high='#000000')

.given_N = .data %>>%
    filter(N == 0) %>>%
    mutate(expected = expected / sum(p))

.given_N %>>% filter(observed > expected)

.range = range(.given_N$expected) %>>% (? .)
.given_N %>>%
    ggplot(aes(expected, observed))+
    geom_point(aes(colour=G+C), alpha=0.4)+
    geom_line(data=data_frame(expected=.range, observed=.range), colour='#FF3300')+
    theme_bw()+
    scale_colour_gradient(low='#00CCCC', high='#000000')

.given_N %>>%
    arrange(observed/expected)

.given_N %>>%
    arrange(desc(observed/expected))

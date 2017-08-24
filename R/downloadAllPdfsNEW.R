library(rvest)
library(xml2)
library(magrittr)
rm(list=ls())

consistentTree  <-  function (...) {
    sptable  <-  readAndParseTree(...)
    while (length(sptable) == 0) {
        sptable  <-  readAndParseTree(...)
    }
    sptable
}

readAndParseTree  <-  function (url) {
    sptree  <-  try(xml2::read_html(url), silent = TRUE)
    while (class(sptree)[1] == 'try-error') {
        sptree  <-  try(xml2::read_html(url), silent = TRUE)
    }
    sptree
}

downloadAllPdfs  <-  function (myDir = getwd(), ...) {
    googleScholarUrls  <-  produceGoogleScholarLinks(...)
    htmlPages          <-  lapply(googleScholarUrls, consistentTree)
    allUrls            <-  unlist(lapply(htmlPages, getPdfUrls))
    sapply(allUrls, downloadPdf, myDir = myDir)
}

produceGoogleScholarLinks  <-  function (searchVector, maxOfPages) { 
    scholarBaseHtml    <-  'https://scholar.google.com.au/scholar?start='
    # GoogleScholar only shows a maximum of 20 per page
    startresults       <-  seq(0, 10 * maxOfPages, 10)
    startresults       <-  startresults[-length(startresults)]
    connectingBit      <-  '&q='
    searchBox          <-  paste0(searchVector, collapse = '+')
    endingBit          <-  '&hl=en&as_sdt=0,5'
    paste0(scholarBaseHtml, startresults, connectingBit, searchBox, endingBit)
}

getPdfUrls  <-  function (htmlPage) {
    urls  <-  htmlPage %>% rvest::html_nodes('a') %>% rvest::html_attr('href')
    urls[grep('*[.]pdf$', urls)]
}

downloadPdf  <-  function (urls, myDir, verbose = TRUE) {
    splits   <-  strsplit(urls, '/')[[1]]
    fileOut  <-  file.path(myDir, splits[length(splits)])
    if (verbose) {
        cat(sprintf('Downloading %s\n', fileOut))
    }
    try(download.file(urls, fileOut), silent = TRUE)
}

dir.create('~/Desktop/tempDownloads')
downloadAllPdfs(myDir = '~/Desktop/tempDownloads', searchVector = c('fish', 'respiration', 'filetype:pdf'), maxOfPages = 1)

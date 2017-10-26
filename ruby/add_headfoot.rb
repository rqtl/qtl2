#!/usr/bin/env ruby
#
# Add header and footer to the Rmd-based vignettes, so they have the
# menubars at the top and the CC-BY at the bottom.

file = ARGV.length > 0 ? ARGV[0] : abort("Give file name as command-line argument")

head_file = "ruby/vignette_header.html"
foot_file = "ruby/vignette_footer.html"

# read header and footer files into single character strings
head = File.readlines(head_file).join
foot = File.readlines(foot_file).join

f = File.open(file)
title = "" # to contain the title
description = "" # to contain the title
author = "" # to contain the title
input = "" # to contain the middle part of the file
passed_header = false
f.each_line do |line|
    if /<title>(.*)<\/title>/ =~ line
        title = $1;
    end

    if /<meta name="description" content="([^"]*)">/ =~ line
        description = $1;
    end

    if /<meta name="author" content="([^"]*)">/ =~ line
        author = $1;
    end

    next if /<h1 class="title toc-ignore">/ =~ line
    next if /<meta/ =~ line

    break if /<\/body>/ =~ line

    input += line if passed_header

    passed_header=true if /<body>/ =~ line
end

abort("Error: Can't find the title") if title == ""

# stick the title into the head string
head.gsub!(/{{TITLE}}/, title)
head.gsub!(/{{DESCRIPTION}}/, description)
head.gsub!(/{{AUTHOR}}/, author)

# close the file being read
f.close

# write over the file with no headers and footers
f = File.open(file, "w")
f.print(head + input + foot)

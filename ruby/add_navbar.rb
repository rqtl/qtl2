#!/usr/bin/env ruby
#
# Add navigation bar to the Rmd-based vignettes, so they have the
# menubars at the top and the CC-BY at the bottom.

file = ARGV.length > 0 ? ARGV[0] : abort("Give file name as command-line argument")

head_file = "ruby/vignette_head.html"
navbar_file = "ruby/vignette_navbar.html"

# read header and footer files into single character strings
head = File.readlines(head_file).join
navbar = File.readlines(navbar_file).join

f = File.open(file)
title = "" # to contain the title
description = "" # to contain the title
author = "" # to contain the title
header = "" # to contain the head part of the file
input = "" # to contain the middle part of the file

passed_head = false

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

    next if /<meta/ =~ line

    if /<\/head>/ =~ line
        passed_head = true
        next
    end

    if /<body>/ =~ line
        line += navbar
    end

    if /<\/body>/ =~ line
        line = "<style>code { background-color: white; }</style>\n" + line;
    end

    if passed_head
        input += line
    else
        header += line
    end

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
f.print(header + head + input)

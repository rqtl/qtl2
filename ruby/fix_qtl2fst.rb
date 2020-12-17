#!/usr/bin/env ruby
#
# Replace the yaml header in qtl2fst.Rmd
# + add CC BY at the bottom

file = ARGV.length > 0 ? ARGV[0] : abort("Give file name as command-line argument")

within_head = false
passed_head = false

contents = ""

f = File.open(file, "r")

f.each_line do |line|
    if passed_head
        if /^### / =~ line
            line.sub!(/^### /, "## ")
        end
        contents += line
    else
        if /^\-\-\-/ =~ line
            if !within_head
                within_head = true
            else
                contents += "output:\n    html_document:\n        toc: true\n        toc_float: true\n"
                passed_head = true
            end
            contents += line
        elsif /title/ =~ line
            contents += line
        end
    end
end

# close the file being read
f.close


contents += "\n\n---\n\n"
contents += "[![CC BY](https://i.creativecommons.org/l/by/3.0/88x31.png)](https://creativecommons.org/licenses/by/3.0/)\n"
contents += "[Karl Broman](https://kbroman.org)\n"


# write over the file with no headers and footers
f = File.open(file, "w")
f.print(contents)

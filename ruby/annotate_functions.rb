#!/usr/bin/env ruby
#
# Annotation all R/qtl2 functions
#
# If pages/rqtl2_functions.md doesn't exist, create it using the S3 methods + exported functions
# in the R/qtl2 NAMESPACE file
#
# If it does exist, check that all of the exported functions are there. Anything missing should
# be appended to the end (in the "Other" section).

mdfile = "pages/rqtl2_functions.md"

# create the pages/annotations file
if !File.exists?(mdfile)
    f = File.open(mdfile, 'w')
    f.print("---\n")
    f.print("layout: page\n")
    f.print("title: R/qtl2 functions\n")
    f.print("description: Annotated/categorized list of functions in R/qtl2\n")
    f.print("---\n\n")
    f.print("### Other functions\n\n")
    f.close()
end

# read the R/qtl2 NAMESPACE file to grab all function names
namespace_file = "../qtl2/NAMESPACE"
functions = []
f = File.open(namespace_file)
f.each_line do |line|
    if line =~ /S3method\(([^)]+)\)/
        the_func = $1
        the_func.sub!(/,/, '.')
        functions.push(the_func)
    elsif line =~ /export\(([^)]+)\)/
        functions.push($1)
    end
end

print("R/qtl2 has #{functions.length()} functions.\n")

# search for function names in the mdfile
f = File.open(mdfile)
f.each_line do |line|
    if line =~ /\- `([^`]+)`/
        the_func = $1
        functions.delete(the_func)
    elsif line =~ /\- \[`([^`]+)`\]/
        the_func = $1
        functions.delete(the_func)
    end

end
f.close()

print("#{functions.length()} functions to add to pages/rqtl2_functions.md\n")

if functions.length() > 0
    f = File.open(mdfile, "a")
    functions.each do |func|
        f.print("- `#{func}`\n")
    end
    f.close()
end

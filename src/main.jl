include("include.jl")

data = createStructures()

fillInput!(data)
fillPreprocessor!(data)
fillOptions!(data)

if data.options.enumerate
    explicityEnumeration!(data)
else
    singleModel!(data)
end
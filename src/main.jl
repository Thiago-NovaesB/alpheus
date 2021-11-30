include("include.jl")

data = createStructures()

fillInput!(data)
fillOptions!(data)
fillPreprocessor!(data)

if data.options.enumerate
    explicityEnumeration!(data)
    grapfOutputEnum!(data)
else
    singleModel!(data)
    grapfOutput!(data)
    textOutput!(data)
end
include("include.jl")

data = createStructures()

fillInput!(data)
fillOptions!(data)
fillPreprocessor!(data)

if data.options.enumerate
    @time explicityEnumeration!(data)
    grapfOutputEnum!(data)
else
    @time singleModel!(data)
    grapfOutput!(data)
    textOutput!(data)
end
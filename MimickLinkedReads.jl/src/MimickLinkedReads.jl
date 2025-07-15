module MimickLinkedReads

using Base.Threads
using BioSequences
using CodecZlib
using Distributions
using FASTX
using Random
using VariantCallFormat
using Term.Progress

include("Structs.jl")
include("Barcodes.jl")
include("Common.jl")
include("FormatFASTQ.jl")
include("Breakpoints.jl")
include("ProcessFASTA.jl")
#include("ProcessFASTQ.jl")
include("ProcessVariants.jl")
include("Main.jl")
# Write your package code here.

end

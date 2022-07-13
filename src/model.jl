abstract type ModelSetup end

"""
    M = BarotropicModel(::Parameters,
                        ::Constants,
                        ::GeoSpectral,
                        ::HorizontalDiffusion)
The BarotropicModel struct holds all other structs that contain precalculated constants,
whether scalars or arrays that do not change throughout model integration. In contrast to
`ShallowWaterModel` or `PrimitiveEquationModel` it does not contain a `Boundaries` struct
as not needed."""
struct Model{NF<:AbstractFloat} <: ModelSetup
    parameters::SystemParameters
    constants::Constants{NF}
end


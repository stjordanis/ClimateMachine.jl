function vars_state_auxiliary(::SoilModelMoisture, FT)
    @vars begin
        z::FT
        h::FT
        κ::FT
    end
end

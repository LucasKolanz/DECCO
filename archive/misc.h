// Generate a random double from -.5lim to .5lim so that numbers are distributed evenly around 0:
double randDouble(double lim)
{
    return lim * ((double)rand() / (double)RAND_MAX - .5);
}

// Returns a vector within the desired radius, resulting in spherical random distribution
double3 randVec(double lim1, double lim2, double lim3)
{
    double3 vec = make_double3(randDouble(lim1), randDouble(lim2), randDouble(lim3));
    double halfLim = lim1 * .5;
    while (length(vec) > halfLim)
    {
        vec = make_double3(randDouble(lim1), randDouble(lim2), randDouble(lim3));
    }
    return vec;
}

// Rounding
std::string rounder(double value, int digits)
{
    return std::to_string(value).substr(0, digits);
}

// Scientific Notation
std::string scientific(double value)
{
    std::stringstream ss;
    ss << value;
    std::string conversion = ss.str();
    return conversion;
}
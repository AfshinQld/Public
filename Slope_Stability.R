#Afshin.Ghahramani@des.qld.gov.au

cat("\014")

# Slope Stability Analysis using the Limit Equilibrium Method
# 
# The Factor of Safety (FoS) is a crucial parameter in slope stability analysis, representing the ratio of resisting forces to driving forces. 
# The following are the key components and equations used in this analysis:
#
# 1. Factor of Safety (FoS):
#    The Factor of Safety is defined as:
#    FoS = (Resisting Forces) / (Driving Forces)
#    A slope is considered stable if the FoS > 1 and unstable if FoS < 1.
#
# 2. Driving Force (F_d):
#    The driving force due to gravity acting on the soil mass is calculated as:
#    F_d = W * sin(θ)
#    Where:
#    W = weight of the soil mass (kN)
#    θ = slope angle (degrees)
#
# 3. Resisting Force (F_r):
#    The resisting force is calculated as:
#    F_r = cA + (W * cos(θ) - uA) * tan(φ)
#    Where:
#    c = cohesion of the soil (kPa)
#    A = area of the potential slip surface (m²)
#    u = pore water pressure (kPa)
#    φ = internal friction angle (degrees)
#
# 4. Calculation of FoS:
#    The Factor of Safety (FoS) is then given by:
#    FoS = F_r / F_d


# Load necessary library
library(dplyr)

# Function to calculate FoS for a given slope
calculate_FoS <- function(weight, slope_angle_degrees, cohesion, friction_angle_degrees, pore_pressure, slip_surface_area) {
  # Convert slope angle and friction angle to radians
  theta_rad <- slope_angle_degrees * pi / 180
  phi_rad <- friction_angle_degrees * pi / 180
  
  # Calculate driving force (F_d)
  driving_force <- weight * sin(theta_rad)
  
  # Calculate resisting force (F_r)
  resisting_force <- cohesion * slip_surface_area + (weight * cos(theta_rad) - pore_pressure * slip_surface_area) * tan(phi_rad)
  
  # Calculate Factor of Safety (FoS)
  FoS <- resisting_force / driving_force
  return(FoS)
}

# Define the data for different materials
materials_data <- data.frame(
  Material_Type = c("Dump", "Waste Rock", "Residual Material"),
  Cohesion_kPa = c(0, 2, 4),                # Typical values for different materials
  Friction_Angle_Degrees = c(25, 35, 20),     # Typical friction angles for different materials
  Pore_Pressure_kPa = c(0, 5, 2),             # Assuming some typical pore pressure conditions
  Unit_Weight_kNm2 = c(15, 20, 18),           # Unit weight of materials in kN/m²
  stringsAsFactors = FALSE
)

# Define the range of slope angles (in degrees) to evaluate
slope_angles <- seq(5, 45, by=5)  # from 5 to 45 degrees, in 5-degree increments

# Define the slip surface area and slope length for calculations
slip_surface_area <- 100  # in square meters, assuming a constant area for simplicity
slope_length <- 100       # in meters, assuming a constant length for simplicity

# Initialize an empty data frame to store results
results <- data.frame()

# Loop through each material type and slope angle
for (i in 1:nrow(materials_data)) {
  for (angle in slope_angles) {
    # Get the properties for the current material
    material <- materials_data[i, ]
    
    # Calculate the weight of the soil mass
    weight_kN <- material$Unit_Weight_kNm2 * slope_length
    
    # Calculate slope percentage
    slope_percentage <- tan(angle * pi / 180) * 100
    
    # Calculate FoS
    FoS <- calculate_FoS(weight_kN, angle, material$Cohesion_kPa, material$Friction_Angle_Degrees, material$Pore_Pressure_kPa, slip_surface_area)
    
    # Store the results
    results <- rbind(results, data.frame(
      Material_Type = material$Material_Type,
      Slope_Angle_Degrees = angle,
      Slope_Percentage = slope_percentage,
      Cohesion_kPa = material$Cohesion_kPa,
      Friction_Angle_Degrees = material$Friction_Angle_Degrees,
      Pore_Pressure_kPa = material$Pore_Pressure_kPa,
      Weight_kN = weight_kN,
      FoS = FoS
    ))
  }
}

# Save the filtered results to a CSV file in C:\temp
write.csv(results, file = "C:/temp/Slope_Stability_1Slope.csv", row.names = FALSE)

# Ensure all columns are displayed
options(tibble.width = Inf) # Set display width to infinite to show all columns
print(results)


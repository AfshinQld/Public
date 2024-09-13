#Afshin.Ghahramani@des.qld.gov.au

cat("\014")

#To make stability data cube - code under progress
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
#
# The following code implements these calculations for various materials classified by the Unified Soil Classification System (USCS).



# Load necessary libraries
library(dplyr)
library(parallel)

# Function to calculate FoS (Factor of Safety) for a given slope
calculate_FoS <- function(weight, slope_angle_degrees, cohesion, friction_angle_degrees, pore_pressure, slip_surface_area) {
  # Convert slope angle and friction angle to radians for calculation
  theta_rad <- slope_angle_degrees * pi / 180
  phi_rad <- friction_angle_degrees * pi / 180
  
  # Driving Force (F_d): Component of the gravitational force acting down the slope
  driving_force <- weight * sin(theta_rad)
  
  # Resisting Force (F_r): Combination of cohesive and frictional forces
  resisting_force <- cohesion * slip_surface_area + 
    (weight * cos(theta_rad) - pore_pressure * slip_surface_area) * tan(phi_rad)
  
  # Calculate Factor of Safety (FoS): Ratio of resisting force to driving force
  FoS <- resisting_force / driving_force
  return(FoS)
}

# Define the data for different materials based on USCS classification
materials_data <- data.frame(
  Material_Type = c("Dump", "Dump", "Dump", "Dump", "Dump", "Dump", 
                    "Waste Rock", "Waste Rock", "Waste Rock", "Waste Rock", 
                    "Residual Material", "Residual Material", "Residual Material", "Residual Material", "Residual Material"),
  USCS_Soil_Type = c("GP", "SP", "SM", "SC", "GM-GL", "SC-CL", 
                     "GW", "SW", "GM", "GC", 
                     "CL", "CH", "MH", "OL", "OH"),
  Description = c("Poorly graded gravel", 
                  "Poorly graded sand", 
                  "Silty sand", 
                  "Clayey sand", 
                  "Silty gravel", 
                  "Clayey sand with many fines", 
                  "Well-graded gravel, fine to coarse gravel", 
                  "Well-graded sand, fine to coarse sand", 
                  "Silty gravel", 
                  "Clayey gravel",
                  "Clay of low plasticity, lean clay", 
                  "Clay of high plasticity, fat clay", 
                  "Silt of high plasticity, elastic silt", 
                  "Organic silt, organic clay", 
                  "Organic clay, organic silt"),
  Cohesion_kPa = c(0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 20, 25, 5, 10, 10), 
  Friction_Angle_Degrees = c(38, 36, 34, 32, 35, 28, 40, 38, 36, 34, 27, 22, 24, 25, 22),
  stringsAsFactors = FALSE
)

# Define the range of vertical heights (in meters) for the slope
vertical_heights <- seq(2, 20, by=2)  # Vertical heights from 2 m to 20 m in 2 m increments

# Define the range of slope angles (in degrees) to evaluate
slope_angles <- seq(5, 30, by=1)  # Slope angles from 5° to 30°

# Define a typical unit weight for the analysis (kN/m²)
unit_weight_kNm2 <- 18

# Define possible pore pressures (kPa) for the analysis
pore_pressures <- seq(0, 10, by=2)  # Example range from 0 to 10 kPa in 2 kPa increments

# Create a grid of all possible combinations of Cohesion and Friction Angle
cohesion_values <- c(0, 5, 10, 15, 20, 25)
friction_angle_values <- c(20, 25, 30, 35, 40)

# Initialize an empty data frame to store results
results <- data.frame()

# Define the number of cores for parallel processing
num_cores <- detectCores() - 1

# Function to perform calculation for each combination of parameters
process_combination <- function(material, height, angle, cohesion, friction_angle, pore_pressure) {
  # Calculate the slope length (m)
  slope_length <- height / sin(angle * pi / 180)
  
  # Calculate the weight of the soil mass (kN)
  weight_kN <- unit_weight_kNm2 * slope_length
  
  # Calculate slope percentage
  slope_percentage <- tan(angle * pi / 180) * 100
  
  # Calculate FoS for the current slope and material properties
  FoS <- calculate_FoS(weight_kN, angle, cohesion, 
                       friction_angle, pore_pressure, slope_length * unit_weight_kNm2)
  
  # Store the results in a data frame, rounding to 2 decimal places, and slope length rounded to integer
  data.frame(
    Material_Type = material$Material_Type,
    USCS_Soil_Type = material$USCS_Soil_Type,
    Description = material$Description,
    Vertical_Height_m = height,
    Slope_Angle_Degrees = round(angle, 2),
    Slope_Length_m = round(slope_length, 0),  # Rounded to nearest integer
    Slope_Percentage = round(slope_percentage, 2),
    Cohesion_kPa = cohesion,
    Friction_Angle_Degrees = friction_angle,
    Pore_Pressure_kPa = pore_pressure,
    Weight_kN = round(weight_kN, 2),
    FoS = round(FoS, 2)
  )
}

# Create a cluster for parallel processing
cl <- makeCluster(num_cores)

# Export necessary objects and functions to the cluster
clusterExport(cl, list("materials_data", "vertical_heights", "slope_angles", "unit_weight_kNm2", "calculate_FoS", "process_combination", "cohesion_values", "friction_angle_values", "pore_pressures"))

# Perform parallel processing
results <- parLapply(cl, 1:nrow(materials_data), function(i) {
  material <- materials_data[i, ]
  do.call(rbind, lapply(vertical_heights, function(height) {
    do.call(rbind, lapply(slope_angles, function(angle) {
      do.call(rbind, lapply(cohesion_values, function(cohesion) {
        do.call(rbind, lapply(friction_angle_values, function(friction_angle) {
          do.call(rbind, lapply(pore_pressures, function(pore_pressure) {
            process_combination(material, height, angle, cohesion, friction_angle, pore_pressure)
          }))
        }))
      }))
    }))
  }))
})

# Stop the cluster
stopCluster(cl)

# Combine the results from all cores
results <- do.call(rbind, results)

# Filter results to ensure FoS is between 1.5 and 2
results <- results %>% filter(FoS >= 1.5 & FoS <= 2)

# Save the filtered results to a CSV file in C:\temp
write.csv(results, file = "C:/temp/Slope_Stability_Results.csv", row.names = FALSE)

# Ensure all columns are displayed without truncation
options(tibble.width = Inf) # Set display width to infinite
print(results)


""" Program to run all the modules in this project sequentially """

include("functions.jl") # import all constants and functions

""" Compute energy production with solar parameters """

println("\nStarting reaction chains with solar parameters.")
Q, E, chain_E = energyprod(T, ρ)
println("Reaction chains complete.\n")

println("     Energy Output\n-----------------------------")
for key in collect(keys(Q))
   @printf("     %-10s | %10.3f MeV\n", key, (Q[key]/MeV))
end

println("\n     Branch Energy Production\n-------------------------------")
for key in collect(keys(chain_E))
   @printf("     %-10s | %10.3f MeV/kgs\n", key, (chain_E[key]/MeV))
end

neutrinoLoss(Q) # run neutrinoLoss function, which prints the losses

print("\nWould you like to perform sanity check? (y/n): ")
ans = lowercase(readline(stdin))
if ans == "y" || ans == "yes"
   sanitycheck(E) # perform sanity check with computed energies
end

print("\nWould you like to plot the energy production as a function of temperature? (y/n): ")

if lowercase(readline(stdin)) == "y"
   print("\nWould you like to save the figure to the current folder? (y/n): ")
   ans = lowercase(readline(stdin))
   if ans == "y" || ans == "yes"
      plot_T(1000, true)
   elseif ans == "n" || ans == "no"
      plot_T(1000)
   end

   gcf()
end

println("\nEnd Program")

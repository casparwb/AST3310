
import fvis3
from os import listdir
import sys

vis = fvis3.FluidVisualiser()

varlist = ["w", "u", "T", "rho", "e", "P"]

animate = input("Would you like to animate? (y/n) ")
if animate == "no" or animate == "n":
    animate=False
elif animate == "yes" or animate == "y":
    animate = True

    folders = listdir("../data")
    print("Please specify from what folder you would like to animate.")
    folder = input(f"Default folder: disturbed_300s.\nPossible choices: {folders}\n")
    if folder == "" or folder == "default":
        folder = "disturbed_300s"
    elif folder not in folders:
        print("Directory not available. Exiting.")
        sys.exit()

    folder = "../data/" + folder
else:
    print("Please give yes or no. Exiting.")
    sys.exit()



while animate:

    var = input("What would you like to animate? ")
    if var not in varlist:
        print(f"Variable not listed. Please give any in {varlist}")
        var = input("What would you like to animate? ")
        if var not in varlist:
            print("I give up.")
            sys.exit()

    save = input("Would you like to save the animation? (y/n) ").lower()
    if save == "no" or save == "n":
        save = False
        videoname=None
    elif save == "yes" or save == "y":
        save = True
        videoname = folder + f"/{var}"
    else:
        print("Please give yes or no. Exiting.")
        sys.exit()


    vis.animate_2D(var, folder=folder, save=save, title=var, video_name=videoname)

    animate = input("Would you like to animate another? (y/n) ")
    if animate == "no" or animate == "n":
        animate = False
    elif animate == "yes" or animate == "y":
        animate = True
    else:
        print("Please give yes or no. Exiting.")
        sys.exit()

def save_everything(foldername):
    titles = ["Vertical Velocity", "Horizontal Velocity",
             "Temperature", "Density", "Energy", "Pressure"]

    datapath = "../data/" + foldername + "/"
    for var, title in zip(varlist, titles):
        savepath = datapath + var

        vis.animate_2D(var, folder=datapath, save=True, \
                       video_name=savepath, title=title)


    n_snapshots = 10
    if "_" in foldername:
        maxtime = float(foldername.split("_")[1][0:-1])
    else:
        maxtime = float(foldername[0:-1])

    snapshot_times = [i/n_snapshots*(maxtime-1) for i in range(11)]

    for (var, title) in zip(varlist,titles):
        vis.animate_2D(var, folder=datapath, \
                       snapshots=snapshot_times, title=title)


ans = input("Would you like to save everything? (y/n) ")
if ans == "no" or ans == "n":
    print("Goodbye.")
    sys.exit()
elif ans == "yes" or ans == "y":
    folders = listdir("../data")
    print("Please choose from which folder.")
    folder = input(f"Possible choices:{folders}\n")
    if folder not in folders:
        print("Directory not available. Exiting.")
        sys.exit()
    save_everything(folder)
else:
    print("Please give yes or no. Exiting.")
    sys.exit()

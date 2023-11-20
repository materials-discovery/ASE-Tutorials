#Exercise 2: Post processing


from ase.io import Trajectory
import pandas as pd
import numpy as np
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt


# Now that we have generated some MD data let's analyse some Cu-Cu distances

# Read in MD trajectory

# In[8]:



traj = Trajectory('X.traj') #Input your MD file here


# get Cu-Cu distances for atoms [0] to [1],[2] from the MD traj and append to list

# In[9]:


bond1 = []
for atoms in traj[0:200]:
    dist1 = atoms.get_distance(0,1)
    bond1.append(dist1)


# In[10]:


bond2 = []
for atoms in traj[0:200]:
    dist2 = atoms.get_distance(0,2)
    bond2.append(dist2)


# ## Prepare data for plotting

# Here we want to open the md logfile as a pandas dataframe

# In[11]:



data = pd.read_csv("md.log", sep='\s{2,}',header=None, nrows=200, skiprows=[0])
data = pd.DataFrame(data)


# create variables with total energy and Cu-Cu distance

# In[12]:


z = data[1].tolist()
x = bond1
y = bond2


# create a dataframe for plotting 3D data

# In[13]:


df = pd.DataFrame(list(zip(x,y,z)), columns=list('XYZ'))
print(df)

# ## Plot a 2D PES

# Now let's try visualising our PES as a 2D contour plot

# In[ ]:


plt.tricontourf(df["X"], df["Y"], df["Z"],levels=10, cmap='plasma')
plt.colorbar()
plt.ylabel('Cu0-Cu1/ Å')
plt.xlabel('Cu0-Cu2/ Å')


# ## Plot density of scatter points

# plot density of scatter points for Cu-Cu distances

# In[ ]:





# In[ ]:


def myplot(x, y, s, bins=1000):
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent
fig, axs = plt.subplots(1, 2, constrained_layout=True)
sigmas = [0, 64]
for ax, s in zip(axs.flatten(), sigmas):
    if s == 0:
        ax.plot(df['X'], df['Y'], 'k.', markersize=5)
        ax.set_title("Scatter plot")
    else:
        img, extent = myplot(df['X'], df['Y'], s)
        ax.imshow(img, aspect="auto", extent=extent, origin='lower', cmap=cm.jet)
        ax.set_title("Density of scatter points")
plt.ylabel('Cu0-Cu1/ Å')
plt.xlabel('Cu0-Cu2/ Å')
plt.show()


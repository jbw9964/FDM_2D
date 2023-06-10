import numpy as np
import matplotlib.pyplot as plt

from matplotlib.cm import coolwarm
from numpy import pi, sin, cos

__all__ = ["Bound_node", "Mesh"]

class Bound_node : 
    pass

class Mesh : 
    def __init__(self, X, Y, dx, dy) :
        self.x_space = np.arange(0, X + dx, dx)
        self.y_space = np.arange(0, Y + dy, dy)
        self.__dx = dx
        self.__dy = dy
        self.__valid_node_num = 0
        self.__valid_node_index = None

        self.xx, self.yy = np.meshgrid(self.x_space, self.y_space, indexing="ij")
        self.grid = np.zeros(shape=(len(self.x_space), len(self.y_space)), dtype=Bound_node)
        self.is_set = False
    
    def return_index(self, x=None, y=None) : 
        if x != None and y != None : 
            x_index = np.isclose(self.x_space, x, atol=self.__dx)
            y_index = np.isclose(self.y_space, y, atol=self.__dy)
            x_index = np.where(x_index)
            y_index = np.where(y_index)
            
            x_return = x_index[0][0]
            y_return = y_index[0][0]
            del x_index, y_index

            return x_return, y_return
        elif x != None : 
            x_index = np.isclose(self.x_space, x, atol=self.__dx)
            x_index = np.where(x_index)

            x_return = x_index[0][0]
            del x_index

            return x_return
        else : 
            y_index = np.isclose(self.y_space, y, atol=self.__dy)
            y_index = np.where(y_index)

            y_return = y_index[0][0]
            del y_index

            return y_return

    def set_rect(self, x_mid_cord, L, H) : 
        x_lower = self.return_index(x_mid_cord - L / 2)
        x_upper = self.return_index(x_mid_cord + L / 2)
        y_lower = self.return_index(x=None, y=0)
        y_upper = self.return_index(x=None, y=0 + H)
        
        self.grid[x_lower:x_upper + 1, y_lower:y_upper + 1] = 1
        self.__valid_node_num = np.count_nonzero(self.grid)
        self.__valid_node_index = np.where(self.grid != 0)
        self.grid[:,:] = 0

        self.grid[x_lower, y_lower:y_upper + 1] = Bound_node()
        self.grid[x_upper, y_lower:y_upper + 1] = Bound_node()
        self.grid[x_lower:x_upper + 1, y_lower] = Bound_node()
        self.grid[x_lower:x_upper + 1, y_upper] = Bound_node()

        self.is_set = True

        del x_lower, x_upper, y_lower, y_upper
    
    def set_half_circle(self, x_mid_cord, R) : 
        # center = (x_mid_cord, 0)
        x_low = self.return_index(x_mid_cord - R)
        x_high = self.return_index(x_mid_cord + R)

        theta = np.linspace(0, pi, 1000)

        x_cord = R * cos(theta) + x_mid_cord
        y_cord = R * sin(theta)

        bound_node_list = []
        bound_node_x = []
        bound_node_y = []
        for x, y in zip(x_cord, y_cord) : 
            x = self.return_index(x=x)
            y = self.return_index(y=y)

            if [x, y] in bound_node_list : 
                continue
            else : 
                bound_node_list.append([x, y])
                bound_node_x.append(x)
                bound_node_y.append(y)
                self.grid[x,:y+1] = 1

        self.__valid_node_num = np.count_nonzero(self.grid)
        self.__valid_node_index = np.where(self.grid != 0)

        for x, y in bound_node_list : 
            self.grid[x-1:x+2, y-1:y+2] = 0
            self.grid[x,:y+1] = 0
        
        self.grid[x_low:x_high + 1, 0] = Bound_node()
        self.grid[tuple([bound_node_x, bound_node_y])] = Bound_node()
        self.is_set = True
        
        del theta, x_cord, y_cord

        del bound_node_list, bound_node_x, bound_node_y
    
    def set_triangular(self, x_mid_cord, base, height) : 

        # left
        # y = ax + b
        # a * (x_low) + b = 0
        # a * (x_mid_index) + b = y_high
        # y_high = a * (x_mid_index - x_low)
        # a = y_high / (x_mid_index - x_low)
        # b = -a * x_low = -1 * (y_high) * (x_low) / (x_low - x_mid_index)

        # right
        # y = ax + b
        # a * (x_high) + b = 0
        # a * (x_mid_index) + b = y_high
        # y_high = a * (x_mid_index - x_high)
        # a = y_high / (x_mid_index - x_high)
        # b = -a * x_high = -1 * (y_high) * (x_high) / (x_high - x_mid_index)

        x_low_index = self.return_index(x_mid_cord - base / 2)
        x_mid_index = self.return_index(x_mid_cord)
        x_high_index = self.return_index(x_mid_cord + base / 2)

        y_high_index = self.return_index(y=height)

        bound_node_list = []
        bound_node_x = []
        bound_node_y = []

        a = y_high_index / (x_mid_index - x_low_index)
        b = -1 * a * x_low_index
        prev_y = 0
        for x in range(x_low_index, x_mid_index + 1) : 
            y = round(a * x + b)
            if [x, y] in bound_node_list : 
                continue
            else : 
                for y_lower in range(prev_y, y + 1) : 
                    bound_node_list.append([x, y_lower])
                    bound_node_x.append(x)
                    bound_node_y.append(y_lower)
                self.grid[x,:y+1] = 1
                prev_y = y

        a = y_high_index / (x_mid_index - x_high_index)
        b = -1 * a * x_high_index
        prev_y = y_high_index
        for x in range(x_mid_index, x_high_index + 1) : 
            y = round(a * x + b)
            if [x, y] in bound_node_list : 
                continue
            else : 
                for y_lower in range(prev_y, y - 1, -1) : 
                    bound_node_list.append([x, y_lower])
                    bound_node_x.append(x)
                    bound_node_y.append(y_lower)
                prev_y = y
                self.grid[x,:y+1] = 1
        
        self.__valid_node_num = np.count_nonzero(self.grid)
        self.__valid_node_index = np.where(self.grid != 0)

        for x, y in bound_node_list : 
            self.grid[x-1:x+2, y-1:y+2] = 0
            self.grid[x,:y+1] = 0
        
        self.grid[x_low_index:x_high_index+1, 0] = Bound_node()
        self.grid[tuple([bound_node_x, bound_node_y])] = Bound_node()
        self.is_set = True

        del bound_node_list, bound_node_x, bound_node_y

    def plot_bound_node(self, elev=45, azim=190, alpha=.9) : 
        bound_index = np.where(self.grid != 0)
        temp_grid = np.zeros_like(self.grid, dtype=float)
        
        index_x = bound_index[0] / len(self.x_space) * self.x_space.max()
        index_y = bound_index[1] / len(self.y_space) * self.y_space.max()
        
        fig, axes = plt.subplots(ncols=1, subplot_kw={"projection":"3d"})
        axes.view_init(elev=elev, azim=azim)
        axes.plot_surface(self.xx, self.yy, temp_grid, alpha=alpha, cmap=coolwarm)
        axes.plot(index_x, index_y, 0, 'o', label="Bound_node", color="#ff7f0e")
        axes.set_xlabel("X")
        axes.set_ylabel("Y")
        axes.legend()
        plt.show()

        del temp_grid, bound_index, fig, axes, index_x, index_y

    def plot_bound_node(self, elev=45, azim=190, alpha=.9) : 
        bound_index = np.where(self.grid != 0)
        temp_grid = np.zeros_like(self.grid, dtype=float)
        
        index_x = bound_index[0] / len(self.x_space) * self.x_space.max()
        index_y = bound_index[1] / len(self.y_space) * self.y_space.max()
        
        fig, axes = plt.subplots(ncols=1, subplot_kw={"projection":"3d"})
        axes.view_init(elev=elev, azim=azim)
        axes.plot_surface(self.xx, self.yy, temp_grid, alpha=alpha, cmap=coolwarm)
        axes.plot(index_x, index_y, 0, 'o', label="Bound_node", color="#ff7f0e")
        axes.set_xlabel("X")
        axes.set_ylabel("Y")
        axes.legend()
        plt.show()

        del temp_grid, bound_index, fig, axes, index_x, index_y

    def plot_2d(self) : 
        bound_index = np.where(self.grid != 0)
        
        index_x = bound_index[0] / len(self.x_space) * self.x_space.max()
        index_y = bound_index[1] / len(self.y_space) * self.y_space.max()
        
        point_x = [0, self.x_space.max(), self.x_space.max(), 0, 0]
        point_y = [0, 0, self.y_space.max(), self.y_space.max(), 0]
        
        plt.plot(point_x, point_y, label="Grid")
        plt.plot(index_x, index_y, 'o', label="Bound_node", color="#ff7f0e")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend()
        plt.show()
        
        del bound_index, index_x, index_y, point_x, point_y

    def get_node_num(self) : 
        if not self.is_set : 
            return 0
        
        return self.__valid_node_num

    def get_valid_node_index(self) : 
        return self.__valid_node_index
    
    def return_grid(self) : 
        return self.grid

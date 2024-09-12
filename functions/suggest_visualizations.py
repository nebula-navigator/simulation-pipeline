
from .plot_histogram import plot_histogram
from .plot_scatter import plot_scatter
from .plot_scatter_3d import plot_scatter_3d
from .plot_boxplot import plot_boxplot
from .plot_pairplot import plot_pairplot
from .plot_violin import plot_violin
from .plot_heatmap import plot_heatmap
from .plot_3d_positions import plot_3d_positions
from .plot_event_frequencies import plot_event_frequencies
from .plot_cumulative_number_vs_stellar_types import plot_cumulative_number_vs_stellar_types
from .plot_position_vs_stellar_type import plot_position_vs_stellar_type
from .plot_cumulative_counts_vs_r import plot_cumulative_counts_vs_r
from .plot_distribution_of_masses import plot_distribution_of_masses
from .plothr import plothr
from .blue_straggler import blue_straggler
import config



def suggest_visualizations(current_data):
    while True:
        data=current_data
        print("\nSuggested visualizations:")
        print("1. Histogram: For visualizing the distribution of a single continuous variable.")
        print("2. Scatter Plot: For visualizing the relationship between two continuous variables, grouped by 2 color axes")
        print("3. 3D Scatter Plot: For visualizing the relationship between three continuous variables.")
        print("4. Box Plot: For visualizing the distribution of a continuous variable, possibly grouped by a categorical variable.")
        print("5. Pair Plot: For visualizing the relationships between multiple continuous variables.")
        print("6. Violin Plot: For visualizing the distribution of a continuous variable, possibly grouped by a categorical variable.")
        print("7. Heatmap: For visualizing the correlation matrix between multiple continuous variables.")
        print("8. 3D Visualization of Positions with Interactive Plot: For visualizing the positions of objects in the cluster in 3D space with interactive features.")
        print("9. Bar Plot of Event Frequencies: For visualizing the frequency of different events (takes single column)")
        print("10. Cumulative Events vs Stellar Types ( bar plot of cumulative sum of each stellar type)")
        print("11. Position vs Stellar Type (radial distribution of each stellar type")
        print("12. Stellar type distribution by mass and radial position. A powerful plot for multiple distribution analysis")
        print("13. Cumulative sum of each stellar type vs radial position line plot (suggested plot from reference thesis)")
        print("14. Plot HR diagram")
        print("15. Plot Blue Straggler cumulative counts vs r")
        print("16. Return to Main Menu")

        viz_choice = input("Enter your choice (1-16): ")

        if viz_choice == '1':
            column_name = input("Enter the column name for histogram: ")
            plot_histogram(data,column_name)
        elif viz_choice == '2':
            x_column = input("Enter the x column name for scatter plot: ")
            y_column = input("Enter the y column name for scatter plot: ")
            color_column = input("Enter the column name for continuous variable (optional): ")
            color_min = input("Enter the minimum value for continuous variable range (optional): ")
            color_max = input("Enter the maximum value for continuous variable range (optional): ")
            categorical_column = input("Enter the column name for categorical variable (optional): ")
            plot_scatter(x_column, y_column,data, color_column if color_column else None, color_min if color_min else None, color_max if color_max else None, categorical_column if categorical_column else None)
        elif viz_choice == '3':
            x_column = input("Enter the x column name for 3D scatter plot: ")
            y_column = input("Enter the y column name for 3D scatter plot: ")
            z_column = input("Enter the z column name for 3D scatter plot: ")
            color_column = input("Enter the column name for color (optional): ")
            plot_scatter_3d(x_column, y_column, z_column,data, color_column if color_column else None)
        elif viz_choice == '4':
            column_name = input("Enter the column name for box plot: ")
            plot_boxplot(column_name)
        elif viz_choice == '5':
            columns_list = input("Enter the column names for pair plot (comma-separated): ").split(',')
            plot_pairplot(data,[col.strip() for col in columns_list])
        elif viz_choice == '6':
            column_name = input("Enter the column name for violin plot: ")
            plot_violin(column_name, data)
        elif viz_choice == '7':
            plot_heatmap(data)
        elif viz_choice == '8':
            color_column = input("Enter the column name for color (optional): ")
            color_min = input("Enter the minimum value for color range (optional): ")
            color_max = input("Enter the maximum value for color range (optional): ")
            plot_3d_positions(data, color_column if color_column else None, color_min if color_min else None, color_max if color_max else None)
        elif viz_choice == '9':
            column_name = input("Enter the column name for bar plot of event frequencies: ")
            plot_event_frequencies(data,column_name)
        elif viz_choice == '10':
            include = input("Include 'single', 'binary', or 'both' stars: ").strip().lower()
            plot_cumulative_number_vs_stellar_types(data,include)
        elif viz_choice == '11':
            plot_position_vs_stellar_type(data)
        elif viz_choice == '12':
            x_column = input("Enter the x column name for distribution of masses: ")
            y_column = input("Enter the y column name for distribution of masses: ")
            plot_distribution_of_masses(data,x_column, y_column)
        elif viz_choice == '13':
            plot_cumulative_counts_vs_r(data)
        elif viz_choice == '14':
            plothr(data)
        elif viz_choice == '15':
            blue_straggler(data)
        elif viz_choice == '16':
            break
        else:
            print("Invalid choice. Please enter a number between 1 and 14.")
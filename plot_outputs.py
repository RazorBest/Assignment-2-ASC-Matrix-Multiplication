import matplotlib.pyplot as plt

DATA_FILE_1 = 'output_graph_neopt'
DATA_FILE_2 = 'output_graph_opt_m'
DATA_FILE_3 = 'output_graph_blas'

def parse_file(filepath):
    x = []
    y = []

    with open(filepath, 'r') as f:
        for line in f.readlines():
            line = line.split('=')
            x.append(float(line[2].split(':')[0]))
            y.append(float(line[3].split('\n')[0]))

    return (x, y)

def main():
    x1, y1 = parse_file(DATA_FILE_1)
    x2, y2 = parse_file(DATA_FILE_2)
    x3, y3 = parse_file(DATA_FILE_3)

    fig, ax = plt.subplots()
    ax.plot(x1, y1, 'o-', linewidth=2.0, label='Times for unoptimised version', 
            color='green')
    plt.annotate(f'N={int(x1[-1])}\nT={y1[-1]:.2f}s', (x1[-1], y1[-1]),
            (x1[-1] - 300, y1[-1] - 10),
            color='green')

    ax.plot(x2, y2, 'o-', linewidth=2.0, label='Times for optimised version',
            color='blue')
    plt.annotate(f'N={int(x2[-1])}\nT={y2[-1]:.2f}s', (x2[-1], y2[-1]),
            (x2[-1] - 100, y2[-1] + 16),
            color='blue')

    ax.plot(x3, y3, 'o-', linewidth=2.0, label='Times for BLAS version',
            color='black')
    plt.annotate(f'N={int(x3[-1])}\nT={y3[-1]:.2f}s', (x3[-1], y3[-1]),
            (x3[-1] - 90, y3[-1] + 8),
            color='black')


    plt.title('Running times for different implementations for B*A*A^t + B^t*B')
    ax.set_ylabel('Time (s)')
    ax.set_xlabel('Size of matrix')
    ax.legend()

    plt.show()


if __name__ == '__main__':
    main()

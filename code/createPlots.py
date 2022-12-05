import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    def plot_nSeqs():
        gibb = np.loadtxt('data/nSeqs_Gibb_Runtime.txt')[0, :]
        meme = np.loadtxt('data/nSeqs_MEME_Runtime.txt')
        x = list(range(10, 51, 2))

        fig = plt.figure()
        plt.plot(x, gibb, 'sb-')
        plt.plot(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Number of Sequences (Linear Scale)')
        plt.xlabel('Number of Sequences')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/nSeqs_reg.jpeg', dpi=500)

        fig = plt.figure()
        plt.semilogy(x, gibb, 'sb-')
        plt.semilogy(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Number of Sequences (LogY Scale)')
        plt.xlabel('Number of Sequences')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/nSeqs_log.jpeg', dpi=500)

    def plot_seqL():
        gibb = np.loadtxt('data/seqL_Gibb_Runtime.txt')[0, :]
        meme = np.loadtxt('data/seqL_MEME_Runtime.txt')
        x = list(range(10, 51, 2))

        fig = plt.figure()
        plt.plot(x, gibb, 'sb-')
        plt.plot(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Sequence Length (Linear Scale)')
        plt.xlabel('Sequence Length')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/seqL_reg.jpeg', dpi=500)
        
        fig = plt.figure()
        plt.semilogy(x, gibb, 'sb-')
        plt.semilogy(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Sequence Length (LogY Scale)')
        plt.xlabel('Sequence Length')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/seqL_log.jpeg', dpi=500)

    def plot_motifL():
        gibb = np.loadtxt('data/motifL_Gibb_Runtime.txt')[0, :]
        meme = np.loadtxt('data/motifL_MEME_Runtime.txt')
        x = list(range(5, 21, 1))

        fig = plt.figure()
        plt.plot(x, gibb, 'sb-')
        plt.plot(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Motif Length (Linear Scale)')
        plt.xlabel('Motif Length')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/motifL_reg.jpeg', dpi=500)

        fig = plt.figure()
        plt.semilogy(x, gibb, 'sb-')
        plt.semilogy(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Motif Length (LogY Scale)')
        plt.xlabel('Motif Length')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/motifL_log.jpeg', dpi=500)
    
    def plot_nMuts():
        gibb = np.loadtxt('data/nMuts_Gibb_Runtime.txt')[0, :]
        meme = np.loadtxt('data/nMuts_MEME_Runtime.txt')
        x = list(range(0, 8, 1))

        fig = plt.figure()
        plt.plot(x, gibb, 'sb-')
        plt.plot(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Number of Mutations (Linear Scale)')
        plt.xlabel('Number of Mutations')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/nMuts_reg.jpeg', dpi=500)

        fig = plt.figure()
        plt.semilogy(x, gibb, 'sb-')
        plt.semilogy(x, meme, 'or-')
        plt.legend(['Gibbs', 'MEME'])
        plt.title('Runtime by Number of Mutations (LogY Scale)')
        plt.xlabel('Number of Mutations')
        plt.ylabel('Runtime (seconds)')
        plt.xticks(x)
        fig.tight_layout()
        fig.savefig('images/nMuts_log.jpeg', dpi=500)

    plot_nSeqs()
    plot_seqL()
    plot_motifL()
    plot_nMuts()
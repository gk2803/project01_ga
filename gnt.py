
import random 

POP_SIZE = 10
random.seed(3)
bits = 50
bounds =[[0,10],[0,20],[0,30]]
Pm = 0.2
Pc=0.7
Er=0.1

def decode( bounds, bits, genes):
        real_chromosome=[]
        for i in range(len(bounds)):
            integer = int(''.join(char for char in genes[i]),2)
            real_value = bounds[i][0] + (integer/(2**bits)) * (bounds[i][1] - bounds[i][0])
            real_chromosome.append(real_value) 
        return real_chromosome

def objective_function(t):
    x = t[0]
    y = t[1]
    z = t[2]

    Objective_max = x**3 + y**3 + z**4 +x*y*z
    
    return Objective_max



class Chromosome:

    def __init__(self, genes:list,prob=0,qprob=0):
        self.genes = genes
        self.prob = prob
        self.qprob = qprob
        self.real_genes = decode(bounds, bits, self.genes)
        self.fitness = objective_function(self.real_genes)
        #self.selected = 0
    
    @classmethod
    def rand(cls):
        chrome=[]
        for _ in range(len(bounds)):
            var=[]
            for _ in range(bits):
                if random.random()>=0.5:
                    var.append("0")
                else:
                    var.append("1")
            chrome.append(var)
        return cls(chrome)

    def get_bin(self):
        s=""
        for gene in self.genes:
            s+=''.join(gene)+" "
        return s

    def __str__(self):
        s = "" 
        t = self.real_genes
        
            
        s+= f" x={t[0]}, y={t[1]}, z={t[2]}, prob={self.qprob}, fitness={self.fitness}"
        return s
        
    def get_genes(self):
        return self.genes




class GeneticAlgorithm:
    def __init__(self,size=POP_SIZE):
        self.size=size
        
        self.population = []
        self.flag = False #if true pop has negative values
        for i in range(size):
            self.population.append(Chromosome.rand())
            if not self.flag and self.population[i].fitness<0:
                self.flag = True
                 

    def misc(self, elitism_ratio = Er):
        self.fitness_sum = 0
        self.qprob = 0
        pop = self.population
        self.elites = []

        #elitism
        sorted_pop = sorted(self.population, key=lambda chromosome: chromosome.fitness, reverse=True)
        for i in range(int(elitism_ratio*self.size)):
            self.elites.append(sorted_pop.pop(i))
        self.population = [chromosome for chromosome in self.population if chromosome  in sorted_pop ]
        print(len(self.population))
        #scale fitness values 
        if self.flag:
            
            min_fitness = sorted_pop[-1].fitness
            for chromosome in pop:
                chromosome.fitness-=min_fitness-10 #fitness can't be zero
                self.fitness_sum += chromosome.fitness
        
        for chromosome in pop:
            self.fitness_sum += chromosome.fitness
        for chromosome in pop:
            chromosome.prob = chromosome.fitness/self.fitness_sum
            self.qprob += chromosome.prob
            chromosome.qprob += self.qprob

        


    def selection(self):
        t=[]

        
        for _ in range(self.size):        
            r = random.random()
            for chromosome in self.population:
                x=chromosome.qprob
                if r <= chromosome.qprob:  
                    t.append(chromosome)
                    break
        
        self.population =t
        
        
    
    
    def crossover(self,crossover_rate:float):
        pop = self.population
        children = list()
        for i in range(int(len(pop)/2)):
            parent1 = pop[2*i-1].get_genes()
            parent2 = pop[2*i].get_genes()
            if random.random()<crossover_rate:
                r = random.random()
                for i in range(1,bits):  
                    if r<i/(bits-1):

                        index = i

                        break
                child1 = [parent1[0][:index] +parent2[0][index:], parent1[1][:index] +parent2[1][index:],parent1[2][:index] +parent2[2][index:] ]
                child2 = [parent2[0][:index] +parent1[0][index:], parent2[1][:index] +parent1[1][index:],parent2[2][:index] +parent1[2][index:] ]
                
                children.append(Chromosome(child1))
                children.append(Chromosome(child2))
            else:
                children.append(Chromosome(parent1))
                children.append(Chromosome(parent2))
        self.population = children
        

    def mutation(self, mutation_rate:float):
        pop = self.population
        offsprings = []
        for chromosome in pop:
            z = chromosome.get_genes()
            if random.random() < mutation_rate:
                i = random.randint(0,2)
                j = random.randint(0,bits-1)
                c = z
                #print("prin")
                #print(c)
                if c[i][j] == '1':#flip
                    c[i][j] = '0' 
                else:
                    c[i][j] ='1'
                offsprings.append(Chromosome(c))
                #print('meta')
                #print(c)
            else:
                offsprings.append(Chromosome(z))
                
        self.population = [*offsprings,*self.elites]

    def best(self):
        best_chrom =self.population[0]
        for chromosome in self.population:
            if chromosome.fitness>best_chrom.fitness:
                best_chrom = chromosome
        return best_chrom

    def __str__(self):
        
        s=""
        for chrome in self.population:
            
            #s+=f"Το χρωμόσωμα {chrome.real_genes} επιλέχθηκε {chrome.selected} φορές ενω αναμενόταν να επιλεχθεί {chrome.prob*self.size} φορές \n"
            s+=f"{chrome.real_genes} fitness ={chrome.fitness} \n"
        return s
        
                    

ga = GeneticAlgorithm(10)


for i in range(1000):
    ga.misc(0.1)
    
    ga.selection()
    ga.crossover(0.5)
    ga.mutation(0.3)
    print("generation: ",i,ga.best())
    #print(ga)
    #for chromosome in ga.elites:
    #    print(chromosome)
print(ga)
#
for chrome in ga.elites:
    print(chrome)

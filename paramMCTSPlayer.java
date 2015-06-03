import java.util.Collections;
import java.util.Random;
import java.util.List;
import java.lang.Math;
import java.util.Stack;
import java.util.Arrays;
import java.lang.Long;

/**
* paramMCTSPlayer - a simple player interface for PokerSquares.
* make placements using Monte Carlo Tree Search with UCT
* algorithm adapted from p9 of "A Survey of Monte Carlo Tree Search Methods"
* @Author: Yawen Chen
* Code adapted from Todd W. Neller
**/

public class ParamMCTSPlayer implements PokerSquaresPlayer{

	private final double EXPLORATION_CONST = 1/(Math.sqrt(2)); // Kcocsis and Szepesvari p8 of the survey 
	private final int PRIM_TIME_FACTOR = 10; 
	private final int SIZE = 5; // number of rows/columns in square grid
	private final int NUM_POS = SIZE * SIZE; // number of positions in square grid
	private final int NUM_CARDS = Card.NUM_CARDS; // number of cards in deck

	private int simulationBound = 1000;

	private double explorationConst = EXPLORATION_CONST;
	private int timeFactor = PRIM_TIME_FACTOR;
	private PokerSquaresPointSystem system; // point system
	private int[] plays = new int[NUM_POS]; // positions of plays so far (index 0 through numPlays - 1) recorded as integers using row-major indices.
	private int numPlays = 0; // number of Cards played into the grid so far
	private Card[][] grid = new Card[SIZE][SIZE]; // grid with Card objects or null (for empty positions)
	private int[][] legalPlayLists = new int[NUM_POS][NUM_POS]; // stores legal play lists indexed by numPlays (depth)
    // (This avoids constant allocation/deallocation of such lists)

	private Card[] simDeck = Card.getAllCards(); // a list of all Cards. As we learn the index of cards in the play deck,
	                                             // we swap each dealt card to its correct index.  Thus, from index numPlays 
	//private List<Card> shuffleDeck = Arrays.asList(simDeck);
	private Random randomCard = new Random(); // pseudorandom number generator for cards
	
	public ParamMCTSPlayer(){
	}

	public ParamMCTSPlayer(double explorationConst, int timeFactor){
		this.explorationConst = explorationConst;
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#init()
	 */
	@Override
	public void init() { 
		// clear grid
		for (int row = 0; row < SIZE; row++)
			for (int col = 0; col < SIZE; col++)
				grid[row][col] = null;
		// reset numPlays
		numPlays = 0;
		// (re)initialize list of play positions (row-major ordering)
		for (int i = 0; i < NUM_POS; i++)
			plays[i] = i;
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#getPlay(Card, long)
	 */
	@Override
	public int[] getPlay(Card card, long millisRemaining) {
		/*
		 *  Player chose the move based on the scores gained from simulation
		 *  For each move, many simulated MCTS are performed 
		 *  For each  play simulation, random undrawn cards are drawn in simulation and the paraMCTS player
		 *   picks a play position that maximizes the score (breaking ties randomly).
		 *   After many such plays, the average score per simulated play is computed.  The play with the highest 
		 *     average score is chosen (breaking ties randomly).   
		 */
		int[] playPos = new int[2];

		if (numPlays == NUM_POS-1) { // last play
			for (int i = 0; i<NUM_POS; i++){
				if (grid[i/SIZE][i % SIZE] == null) {
					playPos[0]=i/SIZE;
                    				playPos[1]=i%SIZE;
                   				grid[i/SIZE][i % SIZE]=card;
				}				
			}			
		}

		else{ // if not the last play
			// copy the play positions (row-major indices) that are empty
			int remainingPlays = NUM_POS - numPlays; // ignores triviality of last play to keep a conservative margin for game completion
			System.arraycopy(plays, numPlays, legalPlayLists[numPlays], 0, remainingPlays);
			

			// match simDeck to actual play event; in this way, all indices forward from the card contain a list of 
			//  undealt Cards in some permutation.
			int cardIndex = numPlays;
			while (!card.equals(simDeck[cardIndex]))
				cardIndex ++;
			simDeck[cardIndex] = simDeck[numPlays];
			simDeck[numPlays] = card;



			if (numPlays ==0){ // don't care where we put the first card
				int randIndex = randomCard.nextInt(NUM_POS); // random number in the range of 0 to 25 
				int play = legalPlayLists[numPlays][randIndex];
				// match plays to the actual event by swtiching
				int index = numPlays;
				while (plays[index] != play)
					index++;
				plays[index] = plays[numPlays];
				plays[numPlays] = play;
				playPos[0] = plays[numPlays] / SIZE; // decode it into row and column
				playPos[1] = plays[numPlays] % SIZE;
			}

			else if (numPlays < NUM_POS -1) { // not the forced last play, nor the first play we don't care about
				// compute average time per move evaluation
				// evaluate time factors 
				long millisPerPlay =  timeEvaluation(remainingPlays, millisRemaining, timeFactor);
				long millisPerMoveEval = millisPerPlay / remainingPlays; // dividing time evenly across moves now considered

				long startTime = System.currentTimeMillis();
	       			long endTime = startTime + millisPerPlay;

				// Push cards after cur cards([numPlays+1: NUM_CARDS]) to the stack
				Stack<Card> stackSimDeck = new Stack<Card>(); 
				for (int j =numPlays+1; j<NUM_CARDS; j++){
					stackSimDeck.push(simDeck[j]);

				}			
				Collections.shuffle(stackSimDeck,randomCard);
				stackSimDeck.push(card); //push the card about to play to the top of the stack

				//Create a root node, passing the states (grid) to the object
				UCT_MCTS root = new UCT_MCTS(numPlays, grid);

				//Constraint budget
				while (System.currentTimeMillis()  < endTime) {
					//run simulation with the current grid status
					root.UCTSearch(stackSimDeck);
				}

				//choose the node with the best score usig UCT
				UCT_MCTS bestChild= root.bestChild(EXPLORATION_CONST);
				playPos= bestChild.returnPlayPosition(card);
			}
			makePlay(card, playPos[0], playPos[1]); // make the chosen play 
		}		
		return playPos; // return the chosen play
	}	

	
	public long timeEvaluation(int remainingPlays, long millisRemaining, int timeFactor){
	// compute average time per move evaluation
	// *** can add factors  to give more time to some moves, ued for simulation. 
	// factor range: 5-100 ( 5 is 20%), suggested: 10
		
		long millisPerPlay = millisRemaining / remainingPlays ; 
		/*
		if (timeFactor >= 5){
			millisPerPlay = millisPerPlay *Long.valueOf(( 1+ ((Math.floor(NUM_POS/4)- Math.abs(remainingPlays-Math.floor(NUM_POS/2)))/Math.floor(NUM_POS/4))/timeFactor));
	        	}
	        	*/
	        	return millisPerPlay; 
	}
	

	public void makePlay(Card card, int row, int col) {
		/*match simDeck to event
		int cardIndex = numPlays;
		while (!card.equals(simDeck[cardIndex]))
			cardIndex++;
		simDeck[cardIndex] = simDeck[numPlays];
		simDeck[numPlays] = card;
		*/
		
		// update plays to reflect chosen play in sequence
		grid[row][col] = card;
		int play = row * SIZE + col;
		int j = 0;
		while (plays[j] != play)
			j++;
		plays[j] = plays[numPlays];
		plays[numPlays] = play;
		
		// increment the number of plays taken
		numPlays++;
	}


	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#setPointSystem(PokerSquaresPointSystem, long)
	 */
	@Override
	public void setPointSystem(PokerSquaresPointSystem system, long millis) {
		this.system = system;
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#getName()
	 */
	@Override
	public String getName() {
		return "ParamMCTSPlayer" ;
	}


	/* 
	*Code based on the algorithm on Page9 of A Survey of Monte Carlo Tree Search Methods
	* The UCT Alogorithm
	*/
	private class UCT_MCTS{

		private Card[][] grids; // current states
		private UCT_MCTS[] children; 
		private UCT_MCTS parent;
		private Stack<Card> cardDeck; //ensure no repitition 
		private int actionsLayer; // layer 0-24
		private int sumOfScore = 0; //rewards
		private int timesVisited = 0; // for UCT
		private int childPosition = 0;
		Boolean expanded = false;

		/* Constructor of node of UCT_MCTS */

		// constructor for root node
		private UCT_MCTS(int numActions, Card[][] grids){
			this.grids = grids;
			this.actionsLayer = numActions;
			this.parent = null;
		}

		// constructor for non-root
		private UCT_MCTS(UCT_MCTS parent){
			this.parent = parent;
			this.actionsLayer = parent.actionsLayer + 1;
		}

		private void UCTSearch(Stack<Card> copySimDeck){
			this.cardDeck = (Stack<Card>) copySimDeck.clone();
			double finalSumScore = 0;
			
			// tree policy
			UCT_MCTS toSelect = this.treePolicy();

			//default policy: running the simulation from toSelect
			for (int i =0; i < simulationBound; i++){		
				finalSumScore += toSelect.defaultPolicy(this.cardDeck);
			}
			
			//Backup
			toSelect.backUp(finalSumScore);

		}

		/*inmplementation of tree policy
		*/
		private UCT_MCTS treePolicy(){
			UCT_MCTS toSelect = this;
	
			while (toSelect.children != null){ // when it is not terminal
				int remainingPlays = NUM_POS - toSelect.actionsLayer;
				if (remainingPlays == 0){ //reach to the terminal
					toSelect.children = null;	//set as terminal
					return toSelect;
				}
				else if (toSelect.expanded == false){ 
					Card curCard =this.cardDeck.pop();
					UCT_MCTS[] children = new UCT_MCTS[remainingPlays];
					for (int i=0; i<remainingPlays; i++){ //expansion if not expanded at all
						toSelect.children[i] = toSelect.expand(curCard);
					}
					toSelect.expanded = true; //set expanded as true	
					int randExplore = randomCard.nextInt(remainingPlays+ 1);
					toSelect = toSelect.children[randExplore];
					return toSelect;
				}
				else{
					toSelect = toSelect.bestChild(EXPLORATION_CONST);
					this.cardDeck.pop(); // in order to match the deck to the node about to expand	
				}													
			}

			return toSelect;		
			}


		/*
		* expand the tree and plays the card 
		*/
		private UCT_MCTS expand(Card curCard){	
			UCT_MCTS newChild = new UCT_MCTS(this);

			Card[][] newGrids = copyGrids(this.grids, SIZE); // construct a new grid for the new node
			// choose random untried actions and assign to the new child
			for (int i =this.childPosition; i< SIZE * SIZE; i++){
				if (this.grids[i/SIZE][i%SIZE] == null){
					newGrids[i/SIZE][i%SIZE] = curCard;
					this.childPosition = i + 1;
					break;
				}
			}
			newChild.grids = newGrids;
			return newChild;
		}


		private Card[][] copyGrids(Card[][] gridsToCopy, int size){
			Card[][] newGrid = new Card[size][size];
			// copy the board from parent
			for (int i =0; i<NUM_POS; i++){
				newGrid[i/size][i%size] = gridsToCopy[i/size][i/size];
			}
			return newGrid;
		}

		private double defaultPolicy(Stack<Card> deckForSim){
			Stack<Card> deckLeft = (Stack<Card>) deckForSim.clone();
			Card[][] simGrids = copyGrids(this.grids, SIZE);

			
			// Fill emptiness on the board randomly
			Collections.shuffle(deckLeft, randomCard);
			for (int i=0; i<NUM_POS; i++){
				if (simGrids[i/SIZE][i%SIZE] == null){
					simGrids[i/SIZE][i%SIZE] = deckLeft.pop();			
				}
			}
			double curScore = system.getScore(simGrids);
			return curScore;
		}
		

		private void backUp(double score){
			UCT_MCTS curNode = this;
			while (curNode != null){
				curNode.timesVisited ++;
				sumOfScore += score;
				curNode = curNode.parent;
			}
		}

		private UCT_MCTS bestChild(double c){
			// return the best child of the node/ root node
			UCT_MCTS bestChild = null;
			double bestScoreSoFar = Double.NEGATIVE_INFINITY;

			double bestScore = Double.MIN_VALUE;
			// select the best score by looping through all children
			for (UCT_MCTS child: this.children){
				// p9 of survey, equation in bestChild
				double curScore = child.sumOfScore/child.timesVisited + c*Math.sqrt(2* Math.log(this.timesVisited)/child.timesVisited);
				if (curScore > bestScoreSoFar){
					bestChild = child;
					bestScoreSoFar = curScore;
				}
			}
			return bestChild;
		}

		private int[] returnPlayPosition(Card card){
			int[] playPos = new int[2];
			for (int i= 0; i<SIZE; i++){
				for (int j=0; j<SIZE; j++){
					if (this.grids[i][j] == card){
						playPos[i] = i;
						playPos[j] = j;
					}
				}
			}
			return playPos;		
		}	
	}	
public static void main(String[] args) {
		PokerSquaresPointSystem system = PokerSquaresPointSystem.getAmeritishPointSystem();
		//PokerSquaresPointSystem system = PokerSquaresPointSystem.getRandomPointSystem();
		System.out.println(system);
		new PokerSquares(new ParamMCTSPlayer(), system).play(); // play a single game
	}	
}
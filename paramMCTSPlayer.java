import java.util.Collections;
import java.util.Random;
import java.util.List;
import java.lang.Math;
import java.util.Stack;

/**
* paramMCTSPlayer - a simple player interface for PokerSquares.
* make placements using Monte Carlo Tree Search with UCT
* algorithm adapted from p9 of "A Survey of Monte Carlo Tree Search Methods"
* @Author: Yawen Chen
* Code adapted from Todd W. Neller
**/

public class ParamMCTSPlayer implements PokerSquaresPlayer{

	private final double EXPLORATION_CONST = 1/(math.sqr(2)); // Kcocsis and Szepesvari p8 of the survey 
	private final int PRIM_TIME_FACTOR; 
	private final int SIZE = 5; // number of rows/columns in square grid
	private final int NUM_POS = SIZE * SIZE; // number of positions in square grid
	private final int NUM_CARDS = Card.NUM_CARDS; // number of cards in deck

	private double explorationConst = EXPLORATION_CONST;
	private int timeFactor = PRIM_TIME_FACTOR;
	private PokerSquaresPointSystem system; // point system
	private int[] plays = new int[NUM_POS]; // positions of plays so far (index 0 through numPlays - 1) recorded as integers using row-major indices.
	private int numPlays = 0; // number of Cards played into the grid so far
	private Card[][] grid = new Card[SIZE][SIZE]; // grid with Card objects or null (for empty positions)
	private int[][] legalPlayLists = new int[NUM_POS][NUM_POS]; // stores legal play lists indexed by numPlays (depth)
    // (This avoids constant allocation/deallocation of such lists during the  selections of MC simulations.)

	private Card[] simDeck = Card.getAllCards(); // a list of all Cards. As we learn the index of cards in the play deck,
	                                             // we swap each dealt card to its correct index.  Thus, from index numPlays 
	private List<Card> shuffleDeck = Arrays.asList(simDeck);
	//private Stack<Card> simDeckStack = new Stack<Card>();
	private Random randomCard = new Random(); // pseudorandom number generator for cards											 // onward, we maintain a list of undealt cards for MC simulation.
	private int[][] legalPlayLists = new int[NUM_POS][NUM_POS]; // stores legal play lists indexed by numPlays (depth)
	// (This avoids constant allocation/deallocation of such lists during the selections of MC simulations.)
	
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
		 * Player chose the move based on the scores gained from simulation
		 *   For each move, many simulated MCTS are performed and the (sometimes
		 *     partially-filled) grid is scored.
		 *   For each  play simulation, random undrawn cards are drawn in simulation and the paraMCTS player
		 *     picks a play position that maximizes the score (breaking ties randomly).
		 *   After many such plays, the average score per simulated play is computed.  The play with the highest 
		 *     average score is chosen (breaking ties randomly).   
		 */

		// copy the play positions (row-major indices) that are empty
		int remainingPlays = NUM_POS - numPlays; // ignores triviality of last play to keep a conservative margin for game completion
		System.arraycopy(plays, numPlays, legalPlayLists[numPlays], 0, remainingPlays);

		// match simDeck to actual play event; in this way, all indices forward from the card contain a list of 
		//   undealt Cards in some permutation.
		int cardIndex = numPlays;
		while (!card.equals(simDeck[cardIndex]))
			cardIndex ++;
		simDeck[cardIndex] = simDeck[numPlays];
		simDeck[numPlays] = card;


		if (numPlays ==0){ // don't care where we put the first card
			int ranIndex = random.nextInt(NUM_POS);
			int play = legalPlayLists[numPlays][ranIndex];
			// match plays to the actual event by swtiching
			int index = numPlays;
			while (plays[index] != play)
				index++;
			plays[index] = plays[numPlays];
			plays[numPlays] = play;
		}

		else if (numPlays < 24) { // not the forced last play, nor the first play we don't care about
			// compute average time per move evaluation
			// evaluate time factors 
			long millisPerPlay =  timeEvaluation(remainingPlays, millisRemaining, 10);
			long millisPerMoveEval = millisPerPlay / remainingPlays; // dividing time evenly across moves now considered

			long startTime = System.currentTimeMillis();
       			long endTime = startTime + millisPerPlay;

			shuffleDeck = Arrays.asList(simDeck[numPlays+1:NUM_CARDS]);
			Collections.shuffle(shuffleDeck,random);
			private Stack<Card> copySimDeck = new Stack<Card>(); 
			for (int i =NUM_CARDS-1; i>=numPlays; i--){
				copySimDeck.push(shuffleDeck[i]);
			}
			copySimDeck.push(simDeck[numPlays]);

			//Create a root node 
			UCT_MCTS root = new UCT_MCTS(numPlays, grid);

			//Constraint budget
			while (System.currentTimeMillis()  < endTime) {
				//run simulation with the current grid status
				root.UCTSearch(copySimDeck);
			}
			//choose the node with the best score usig UCT
			root.bestChild(0);

			//clear out
			// clear nodes and grids, code goes here:


				

			/*delete
			double maxAverageScore = Double.NEGATIVE_INFINITY; // maximum average score found for moves so far
			ArrayList<Integer> bestPlays = new ArrayList<Integer>(); // all plays yielding the maximum average score 

			for (int i = 0; i < remainingPlays; i++) { // for each legal play position
				int play = legalPlayLists[numPlays][i];
				long startTime = System.currentTimeMillis();
				long endTime = startTime + millisPerMoveEval; // compute when MC simulations should end
				makePlay(card, play / SIZE, play % SIZE);  // play the card at the empty position
				int simCount = 0;
				int scoreTotal = 0;
				while (System.currentTimeMillis() < endTime) { // perform as many MCTS simulations as possible through the allotted time
					// Perform a Monte Carlo simulation of greedy play to the depth limit or game end, whichever comes first.
					scoreTotal += MCTSPlay();  // accumulate MCTS simulation scores
					simCount++; // increment count of MCTS simulations
				}
				undoPlay(); // undo the play under evaluation
				// update (if necessary) the maximum average score and the list of best plays
				double averageScore = (double) scoreTotal / simCount;
				if (averageScore >= maxAverageScore) {
					if (averageScore > maxAverageScore)
						bestPlays.clear();
					bestPlays.add(play);
					maxAverageScore = averageScore;
				}
			}
			int bestPlay = bestPlays.get(random.nextInt(bestPlays.size())); // choose a best play (breaking ties randomly)
			// update our list of plays, recording the chosen play in its sequential position; all onward from numPlays are empty positions
			int bestPlayIndex = numPlays;
			while (plays[bestPlayIndex] != bestPlay)
				bestPlayIndex++;
			plays[bestPlayIndex] = plays[numPlays];
			plays[numPlays] = bestPlay;
		}
		*/
		int[] playPos = {plays[numPlays] / SIZE, plays[numPlays] % SIZE}; // decode it into row and column
		makePlay(card, playPos[0], playPos[1]); // make the chosen play (not undoing this time)
		return playPos; // return the chosen play
	}

	private long timeEvaluation(int remainingPlays, long millisRemaining, int factor){
	// compute average time per move evaluation
	// *** can add factors  to give more time to some moves, ued for simulation. 
	// factor range: 5-100 ( 5 is 20%), suggested: 10
		long millisPerPlay = millisRemaining / remainingPlays ; 
		if (factor >= 5){
			millisPerPlay = millisPerPlay *( 1+ ((Math.floor(NUM_POS/4)- Math.abs(remainingPlays-Math.floor(NUM_POS/2)))/Math.floor(NUM_POS/4))/factor);
	        	}
	        	return millisPerPlay 
	}

	/**
	 * From the chosen play, perform simulated Card draws with MCTS
	 * and return the resulting grid score.
	 * other prarameters to add: factors for time allocation 
	 * @return resulting grid score after MCTS simulation to given depthLimit
	 */
	private int MCTSPlay2() {
	// ignores this now, code moved to UCT_MCTS class
		return;
	}

	public void makePlay(Card card, int row, int col) {
		// match simDeck to event
		int cardIndex = numPlays;
		while (!card.equals(simDeck[cardIndex]))
			cardIndex++;
		simDeck[cardIndex] = simDeck[numPlays];
		simDeck[numPlays] = card;
		
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

	public void undoPlay() { // undo the previous play
		numPlays--;
		int play = plays[numPlays];
		grid[play / SIZE][play % SIZE] = null;	
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

		private UCTSearch(Stack<Card> copySimDeck){
			this.cardDeck = (Stack<Card>) copySimDeck.clone();
			double finalSumScore;
			
			// tree policy
			UCT_MCTS toSelect = this.treePolicy();

			//default policy: running the simulation from toSelect
			for (i =0; i < simulationBound; i++){		
				finalSumScore += toSelect.defaultPolicy(this.cardDeck);
			}
			
			//Backup
			toSelect.backUp(curScore);

		}

		/*inmplementation of tree policy
		*/
		private UCT_MCTS treePolicy(){
			UCT_MCTS toSelect = this;
	
			while (toSelect.children != null){ // when it is not terminal
				int remainingPlays = NUM_POS - toSelect.actionsLayer;
				if (remainingPlays == 0){ //reach to the terminal
					toSelect.children = null;	//set as terminal
					return;
				}
				else if (toSelect.expanded == false){ 
					Card curCard =this.cardDeck.pop();
					UCT_MCTS[] children = new UCT_MCTS[remainingPlays];
					for (int i; i<remainingPlays; i++){ //expansion if not expanded at all
						toSelect.children[i] = toSelect.expand(curCard);
					}
					toSelect.expanded = true; //set expanded as true	
					int randExplore = math.random() * ((remainingPlays) + 1);
					toSelect = toSelect.children[randExplore];
					return toSelect;
				}
				else{
					toSelect = toSelect.bestChild(EXPLORATION_CONST)
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
			for (i =this.childPosition; i< SIZE * SIZE; i++){
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
			for (i =0; i<NUM_POS; i++){
				newGrid[i/size][i%size] = gridsToCopy[i/size][i/size];
			}
			return newGrid
		}

		private double defaultPolicy(Stack<Card> deckForSim){
			Stack<Card> deckLeft = (Stack<Card>) deckForSim.clone();
			Card[][] simGrids = copyGrids(this.grids, SIZE)
			
			// Fill emptiness on the board randomly
			Collections.shuffle(deckLeft, random);
			for (i=0; i<NUM_POS; i++){
				if (simGrids[i/SIZE][i%SIZE] == null){
					simGrids[i/SIZE][i%SIZE] = deckLeft.pop();			
				}
			}
			double curScore = system.getScore(simGrids)
			return curScore;
		}
		

		private backUp(double score){
			UCT_MCTS curNode = this;
			while (curNode != null){
				curNode.timesVisited ++;
				sumOfScore += score;
				curNode = curNode.parent;
			}
		}

		private UCT_MCTS bestChild(int c){
			// return the best child of the node/ root node
			UCT_MCTS bestChild;
			double bestScoreSoFar = Double.NEGATIVE_INFINITY;

			double bestScore = Double.MIN_VALUE;
			// select the best score by looping through all children
			for (UCT_MCTS child: this.children){
				// p9 of survey, equation in bestChild
				double curScore = child.sumOfScore/child.timesVisited + c*Math.sqrt(2* Math.log(this.timesVisited)/child.timesVisited)
				if (curScore > bestScoreSoFar){
					bestChild = child;
					bestScoreSoFar = curScore
				}
			}
			return bestChild
	}
public static void main(String[] args) {
		PokerSquaresPointSystem system = PokerSquaresPointSystem.getAmeritishPointSystem();
		//PokerSquaresPointSystem system = PokerSquaresPointSystem.getRandomPointSystem();
		System.out.println(system);
		new PokerSquares(new ParamMCTSPlayer(), system).play(); // play a single game
	}	
}s
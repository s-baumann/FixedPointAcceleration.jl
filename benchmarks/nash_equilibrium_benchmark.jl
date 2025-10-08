using BenchmarkTools
using FixedPointAccelerationNext
using LinearAlgebra
using Random

"""
Computing Nash Equilibrium in an N-player game using fixed point iteration.
Each player's strategy is a probability distribution over actions.
"""

struct NPlayerGame{T<:Real}
    n_players::Int
    n_actions::Vector{Int}  # actions per player
    payoffs::Vector{Array{T}}  # payoff tensors for each player
end

function random_game(n_players::Int, actions_per_player::Int; seed=42)
    Random.seed!(seed)
    n_actions = fill(actions_per_player, n_players)

    # Create random payoff tensors
    payoffs = []
    for i in 1:n_players
        # Payoff tensor has dimensions (actions_per_player)^n_players
        dims = tuple(n_actions...)
        payoff_tensor = randn(dims...)
        push!(payoffs, payoff_tensor)
    end

    return NPlayerGame{Float64}(n_players, n_actions, payoffs)
end

function best_response_operator(game::NPlayerGame)
    n_players = game.n_players
    n_actions = game.n_actions
    total_vars = sum(n_actions)

    return function (strategies_flat::Vector{T}) where {T}
        # Reconstruct strategy profiles from flat vector
        strategies = []
        idx = 1
        for i in 1:n_players
            n_act = n_actions[i]
            player_strategy = strategies_flat[idx:(idx + n_act - 1)]
            # Ensure probability simplex (normalize)
            player_strategy = max.(player_strategy, 1e-10)  # Avoid zeros
            player_strategy = player_strategy / sum(player_strategy)
            push!(strategies, player_strategy)
            idx += n_act
        end

        new_strategies_flat = similar(strategies_flat)
        idx = 1

        for player in 1:n_players
            n_act = n_actions[player]
            payoff_tensor = game.payoffs[player]

            # Calculate expected payoffs for each action
            expected_payoffs = zeros(T, n_act)

            for action in 1:n_act
                # Create index tuple with player's action fixed
                expected_payoff = 0.0

                # Iterate over all opponent strategy combinations
                opponent_indices = [1:n_actions[j] for j in 1:n_players if j != player]

                if length(opponent_indices) > 0
                    for opponent_combo in Iterators.product(opponent_indices...)
                        # Build full index tuple
                        full_index = Vector{Int}(undef, n_players)
                        full_index[player] = action

                        opp_idx = 1
                        for j in 1:n_players
                            if j != player
                                full_index[j] = opponent_combo[opp_idx]
                                opp_idx += 1
                            end
                        end

                        # Calculate probability of this outcome
                        prob = 1.0
                        opp_idx = 1
                        for j in 1:n_players
                            if j != player
                                prob *= strategies[j][opponent_combo[opp_idx]]
                                opp_idx += 1
                            end
                        end

                        # Add to expected payoff
                        expected_payoff += prob * payoff_tensor[full_index...]
                    end
                else
                    # Single player case
                    expected_payoff = payoff_tensor[action]
                end

                expected_payoffs[action] = expected_payoff
            end

            # Best response: softmax to ensure smooth iteration
            β = 10.0  # Temperature parameter
            exp_payoffs = exp.(β * expected_payoffs)
            best_response = exp_payoffs / sum(exp_payoffs)

            # Store in flat vector
            new_strategies_flat[idx:(idx + n_act - 1)] = best_response
            idx += n_act
        end

        return new_strategies_flat
    end
end

function create_nash_benchmark(n_players::Int=3, actions_per_player::Int=4)
    game = random_game(n_players, actions_per_player)
    T = best_response_operator(game)

    # Initial strategy: uniform random
    total_vars = sum(game.n_actions)
    x0 = rand(total_vars)

    # Normalize each player's strategy to probability simplex
    idx = 1
    for i in 1:n_players
        n_act = game.n_actions[i]
        x0[idx:(idx + n_act - 1)] =
            x0[idx:(idx + n_act - 1)] / sum(x0[idx:(idx + n_act - 1)])
        idx += n_act
    end

    cfg = FixedPointConfig(; threshold=1e-8, max_iters=1000)

    return T, x0, cfg, game
end

function nash_equilibrium!(suite)
    T, x0, cfg, game = create_nash_benchmark(3, 3)  # 3 players, 3 actions each

    methods = [
        ("Simple", Simple()),
        ("Anderson", Anderson(; m=6)),
        ("Aitken", Aitken()),
        ("MPE", MPE(; period=4)),
        ("RRE", RRE(; period=4)),
        ("VEA", VEA(; period=6)),
        ("SEA", SEA(; period=6)),
    ]

    for (name, method) in methods
        suite["Real-life problem"]["Nash equilibrium"][name] = @benchmarkable solve(
            $T, x0_copy; method=($method), cfg=($cfg)
        ) setup=(x0_copy = copy($x0))
    end
end
